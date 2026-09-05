(function () {
  'use strict';

// App shell. Loads an assembly FASTA (streamed through fasta-index.js in a
// Worker, with marker-gene search) and any number of contig->bin tables
// alongside it — content-sniffed from whatever else was selected, not by
// filename (the tool label per table IS taken from its filename, since
// content alone can't name which binning tool produced it). Two or more
// tables get cross-tool reconciliation: bins matched across tools by
// contig overlap, core/disputed contig sets.
//
// The app is deliberately centred on one workflow: pick a putative MAG
// (the picker table in renderReconciliationCard) to see it, and any MAG
// its disputed contigs also touch, as a contig network
// (buildMagNeighborhood/initMagNetwork); click a contested contig to see
// per-candidate-MAG evidence — GC/coverage vs. that MAG's core contigs,
// marker-gene unique-vs-duplicate status — and decide where it belongs or
// exclude it (renderContigEvidence). Every decision writes into
// `workingAssignment`, the one live session state the MAG picker table,
// the network's leaf colouring, and Export all read from. Earlier
// standalone views (per-tool bin summaries, a scatter-based manual
// reassignment section, whole-assembly QC/redundancy comparisons) were
// removed in that redesign — see docs/HANDOVER.md — in favour of this
// single contig-resolution flow.

const THEME_KEY = 'clann-mag-explorer-theme';

// Session-level interactive-reassignment state — deliberately module-level
// (not passed around as parameters) since it needs to survive across
// re-renders and outlives any single render() call. Reset on every fresh
// assembly load (loadAssembly) — a new session starts clean.
let workingAssignment = new Map();
let workingAssignmentInitialized = false; // guards deriveInitialAssignment from re-running (and losing manual edits) when a filter change rebuilds the DOM

// Which putative MAG's neighborhood the contig network is currently scoped
// to, and which contig within it has an evidence panel open — both null
// until the student picks one from the MAG picker table. Reset on every
// fresh assembly load, same as workingAssignment above.
let selectedMagId = null;
let selectedContigId = null;

// Network scope toggles, both default off: the network starts scoped to
// just the selected MAG's *contended* contigs (the thing worth resolving)
// plus the bare hubs of whatever else contends with them, not every
// uncontended contig the primary MAG has and not every one of a secondary
// MAG's own contigs — pulling in a secondary MAG's full contig set by
// default turned out to be the wrong call (it buries the actual dispute
// under everything else that MAG happens to contain), so both are now
// explicit, off-by-default choices instead. Reset on every fresh load.
let showUncontendedContigs = false;
let showConnectedMagContigs = false;

// Contigs the student has explicitly decided via the evidence panel's
// "Assign here"/"Exclude" buttons — deliberately NOT the same thing as
// "workingAssignment has an entry for this contig", since
// deriveInitialAssignment already seeds workingAssignment with a
// majority-vote default for every voted contig (including disputed ones)
// as a starting point. Without this separate set, every disputed contig
// would render as 'resolved' from the moment of load, before the student
// touched anything — this set is what actually distinguishes "still
// needs a look" from "you decided this."
let decidedContigIds = new Set();

// Sentinel working-assignment "bin" for a contig the student has decided
// belongs to neither of its contended MAGs — a real value in the same
// contigId->binId Map as every other assignment (working-assignment.js's
// reassignContigs doesn't need to know this is special), so it shows up
// like any other bin in the export dropdown/CSV rather than needing its
// own separate tracking structure.
const EXCLUDED_BIN_ID = 'excluded';

// Left-pane filters (brief §Left pane: filters and search) — `currentFilters`
// is the live filter state (src/model/filters.js's shape), and `latest` is
// everything about the current load that filtering must NOT recompute:
// the full record set, the loaded bin tables, the (filter-independent)
// cross-tool reconciliation result and outlier flags, and the lookup
// indices filters.js needs. Filtering only changes which subset of that
// fixed data gets rendered (renderFilteredExplorer), never recomputes
// reconciliation/outliers/marker search themselves — those describe the
// whole loaded dataset, not "the dataset restricted to what's currently
// visible". Reset on every fresh assembly load, same as workingAssignment.
let currentFilters = null;
let latest = null;

// MAG-level filtering (src/model/mag-filters.js) — a second, independent
// filter axis over the cross-tool reconciliation table's own columns
// (Putative MAG, Core, Disputed, Completeness, Redundancy, Tier, supporting
// tool), rather than per-contig properties. This is the left pane's primary
// job in the redesign: narrowing the MAG picker table down to the MAG
// worth selecting next.
let currentMagFilters = null;

// Global thresholds/parameters a student can adjust to see how sensitive
// the tool's derived calls are to where these lines are drawn: MIMAG
// quality-tier cutoffs (bin-summary.js), the cross-tool bin-matching
// overlap threshold (bin-reconciliation.js's minJaccard), and the outlier
// flagging thresholds (outliers.js's computeFlagCount). Everything
// here is cheap, pure-JS recomputation over already-parsed data — no
// worker/marker-search re-run — except minJaccard, which changes bin
// matching itself and so triggers a full recompute (recomputeLatest) at
// the same cost as a fresh load's compute-once phase; every other
// parameter only changes what gets *displayed* from data already in
// `latest`, so those just re-render.
function defaultGlobalParams() {
  return {
    mimag: { ...window.ClannMAG.binSummary.DEFAULT_MIMAG_THRESHOLDS },
    minJaccard: 0.1,
    outlier: { ...window.ClannMAG.outliers.DEFAULT_OUTLIER_PARAMS },
    // Recall-adjustment for completeness/redundancy (docs/scg-blast-
    // verification.md) — see bin-summary.js's DEFAULT_ESTIMATED_RECALL
    // for why this correction exists at all: the built-in marker search
    // is a fast, approximate heuristic (not a profile-HMM search), and an
    // independent BLAST verification measured it only recovers ~74% of
    // genuinely-present marker genes even after this app's own threshold
    // tuning. Reading "families found / 40" as completeness would
    // silently understate every genuinely-complete genome by roughly
    // that same gap; dividing by the measured recall rate corrects for
    // it instead of reporting a number known to be biased low.
    recallRate: window.ClannMAG.binSummary.DEFAULT_ESTIMATED_RECALL,
  };
}
let currentParams = null;

// The complete set of loaded contig->bin tables, exactly as parsed — never
// filtered — alongside which of those tools are currently "active"
// (included in reconciliation). `latest.binTablesByTool` (built by
// recomputeLatest) is always the *active subset* of this; keeping the full
// map separately is what lets a student re-enable an excluded tool without
// re-uploading it. Both reset on every fresh assembly load, `activeTools`
// to every loaded tool (nothing excluded by default).
let allBinTablesByTool = null;
let activeTools = new Set();

/** The active-only view of allBinTablesByTool that recomputeLatest should run on. */
function filterActiveBinTables() {
  if (!allBinTablesByTool) return null;
  const filtered = new Map();
  for (const [tool, assignments] of allBinTablesByTool) {
    if (activeTools.has(tool)) filtered.set(tool, assignments);
  }
  return filtered.size ? filtered : null;
}

// Persists the reconciliation network's chosen layout algorithm across
// filter-triggered re-renders (which rebuild the network's DOM each time —
// see renderFilteredExplorer), same reasoning as workingAssignmentInitialized
// above. Reset on every fresh assembly load.
let networkAlgorithm = 'ring';

// Phase 9 export: the loaded assembly's own File (re-sliced via
// Blob.slice for per-MAG FASTA extraction, never re-uploaded) and its
// parsed records, kept at module scope alongside workingAssignment since
// export reads from the same live session state. Reset on every fresh
// load, same as workingAssignment above — a reload loses the File
// reference regardless (brief's accepted, documented limitation; see
// docs/phase1-investigation.md Phase 1), so there's nothing to persist.
let currentAssemblyFile = null;
let currentRecords = [];

// Broader than the brief's original design (docs/clann-mag-explorer-brief.md
// §Unsaved-work protection, which armed this only after the first contig
// reassignment, deliberately not before, "so it doesn't nag a student
// who's only exploring"). Changed on direct request: any successfully
// loaded assembly now arms the guard, not just an edited one — losing the
// in-memory File reference (needed for per-MAG FASTA extraction, see
// currentAssemblyFile above) on an accidental navigation is disruptive
// even before any reassignment has been made. Gated on currentRecords
// (set once renderContigTable runs after a successful parse — see below),
// not currentAssemblyFile (set as soon as a load is attempted, before
// parsing finishes), so a failed/in-progress load doesn't arm this.
window.addEventListener('beforeunload', (e) => {
  if (currentRecords.length === 0) return;
  e.preventDefault();
  e.returnValue = ''; // required by some browsers to trigger the native confirmation dialog
});

function initTheme() {
  const saved = localStorage.getItem(THEME_KEY);
  if (saved === 'light' || saved === 'dark') {
    document.documentElement.setAttribute('data-theme', saved);
  }
  const btn = document.getElementById('themeBtn');
  btn.addEventListener('click', () => {
    const current = document.documentElement.getAttribute('data-theme');
    const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
    const isDark = current ? current === 'dark' : prefersDark;
    const next = isDark ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    localStorage.setItem(THEME_KEY, next);
  });
}

function showError(message) {
  const err = document.getElementById('err');
  if (!message) { err.style.display = 'none'; err.textContent = ''; return; }
  err.style.display = 'block';
  err.textContent = message;
}

// Header spinner shown whenever any long-running process (currently:
// loadAssembly's parse + marker-gene search) is in flight. A counter, not
// a boolean, so overlapping callers (there's only one today, but this is
// the cheap way to make that not a landmine later) can each set/clear
// independently without one finishing early hiding it for the other.
let busyCount = 0;
function setBusy(active) {
  busyCount = Math.max(0, busyCount + (active ? 1 : -1));
  document.getElementById('hSpinner').style.display = busyCount > 0 ? 'inline-block' : 'none';
}

/** Cheap content-based FASTA check: gzip-aware peek at the first non-blank byte. */
async function looksLikeFasta(file) {
  const { isGzipped } = window.ClannMAG.fastaIndex;
  if (await isGzipped(file)) return true; // trust streamFasta's own decompression + parse to confirm
  const head = await file.slice(0, 1).text();
  return head === '>';
}

function computeN50(lengthsDesc, totalLength) {
  let cumulative = 0;
  for (const len of lengthsDesc) {
    cumulative += len;
    if (cumulative >= totalLength / 2) return len;
  }
  return 0;
}

function formatMarkerHits(markerHits) {
  if (!markerHits) return '<span class="hint" title="Marker-gene assets did not load">n/a</span>';
  if (markerHits.length === 0) return '';
  return markerHits
    .map((h) => `<span title="${h.representativeCount} representatives agreed, score ${h.bestScore}, provenance taxID ${h.provenanceTaxId}">${h.family}</span>`)
    .join(', ');
}

function formatMimagTier(tier) {
  const label = { high: 'High', medium: 'Medium', low: 'Low' }[tier] || tier;
  return `<span class="mimag-tier mimag-tier-${tier}">${label}</span>`;
}

/**
 * Builds the MAG picker's per-MAG numbers as plain data (not HTML) — this
 * is both what renderReconciliationCard renders AND exactly the shape
 * src/model/mag-filters.js's applyMagFilters expects. Core/disputed counts
 * are the original, static cross-tool vote counts (stable across a
 * session, useful for "which MAGs had the most disputes worth resolving").
 * Completeness/redundancy/tier, though, are read off the *live* working
 * assignment (working-assignment.js) rather than the static reconciled
 * core set — so once the student starts resolving disputed contigs via the
 * evidence panel, this table's numbers move with those decisions, the
 * same "every bin-level summary statistic recalculates live" property the
 * old Phase 7 scatter view had, now folded into this one table instead of
 * a separate section. Computed once per render pass in renderFilteredExplorer
 * and threaded through, so filtering doesn't redo this work per keystroke.
 */
function computeMagSummaryData(records, result) {
  const { computeBinSummaries } = window.ClannMAG.binSummary;
  const { assignmentToRows } = window.ClannMAG.workingAssignment;
  const magIds = new Set(result.putativeMags.map((m) => m.magId));
  const liveRows = assignmentToRows(workingAssignment).filter((r) => magIds.has(r.binId));
  const { summaries } = computeBinSummaries(records, liveRows, { thresholds: currentParams.mimag, recallRate: currentParams.recallRate });
  const liveByMagId = new Map(summaries.map((s) => [s.binId, s]));

  return result.putativeMags.map((mag) => {
    const live = liveByMagId.get(mag.magId);
    return {
      magId: mag.magId,
      coreCount: mag.coreContigIds.length,
      disputedCount: mag.disputedContigIds.length,
      liveContigCount: live ? live.contigCount : 0,
      completeness: live ? live.completeness : 0,
      redundancy: live ? live.redundancy : 0,
      tier: live ? live.mimagTier : 'low',
      tools: mag.members.map((m) => m.tool),
      mag,
    };
  });
}

/**
 * Renders the MAG picker table — the sole entry point into the contig
 * network below. Selecting a row (via .mag-picker-select, wired as a
 * delegated document click listener so it survives this card's wholesale
 * innerHTML rebuilds) sets `selectedMagId` and re-renders, which scopes
 * the network to that MAG's neighborhood (see buildMagNeighborhood).
 * `filteredMagIds` (from the left-pane MAG filters — the left pane's job
 * now, per the redesign, is finding which MAG to select here) restricts
 * which rows show; core/disputed are the static cross-tool vote counts,
 * completeness/redundancy/tier are live off the working assignment (see
 * computeMagSummaryData).
 */
function renderReconciliationCard(records, result, magSummaryData, filteredMagIds) {
  const tools = result.tools;

  const magRows = magSummaryData
    .filter((m) => filteredMagIds.has(m.magId))
    .map((m) => {
      const toolCells = tools
        .map((tool) => {
          const member = m.mag.members.find((mem) => mem.tool === tool);
          return member
            ? `<td>${member.binId} <span class="hint">(${member.contigCount.toLocaleString()})</span></td>`
            : '<td class="hint">—</td>';
        })
        .join('');
      const isSelected = m.magId === selectedMagId;
      return `<tr class="mag-picker-row${isSelected ? ' mag-picker-row-selected' : ''}">
        <td><button class="act mag-picker-select" type="button" data-mag-id="${m.magId}">${isSelected ? '● ' : ''}${m.magId}</button></td>
        <td class="num">${m.liveContigCount.toLocaleString()}</td>
        <td class="num">${m.coreCount.toLocaleString()}</td>
        <td class="num">${m.disputedCount.toLocaleString()}</td>
        <td class="num">${m.completeness.toFixed(1)}%</td>
        <td class="num">${m.redundancy.toFixed(1)}%</td>
        <td>${formatMimagTier(m.tier)}</td>
        ${toolCells}
      </tr>`;
    })
    .join('');

  return `
    <div class="card">
      <h3>Cross-tool reconciliation</h3>
      <div class="row-count">${tools.length} tools loaded (${tools.join(', ')}) &middot; ${magSummaryData.length.toLocaleString()} putative MAGs matched by contig overlap (reciprocal best hit, min Jaccard ${currentParams.minJaccard}) &middot; ${filteredMagIds.size.toLocaleString()} of ${magSummaryData.length.toLocaleString()} match the current MAG filters &middot; select a MAG to explore it below. Contigs/Completeness/Redundancy reflect your current working decisions (see Export); Core/Disputed are the original cross-tool vote counts.</div>
      <div class="table-wrap scroll-panel">
        <table class="data-table">
          <thead><tr>
            <th>Putative MAG</th>
            <th title="Contigs currently assigned here in your working decisions">Contigs</th>
            <th title="Contigs every voting tool originally agreed belong to this MAG">Core</th>
            <th title="Contigs originally assigned here by some but not all voting tools">Disputed</th>
            <th title="Recall-adjusted, from the current working assignment — see Thresholds & parameters">Completeness</th>
            <th title="Recall-adjusted, from the current working assignment — see Thresholds & parameters">Redundancy</th>
            <th>Tier</th>
            ${tools.map((t) => `<th>${t}</th>`).join('')}
          </tr></thead>
          <tbody>${magRows}</tbody>
        </table>
      </div>
      <h4>Contig network</h4>
      <div class="row"><label title="Off by default: only contended contigs (voted differently by at least one tool) show. Turn on to also show this MAG's uncontended contigs.">
        <input type="checkbox" id="toggleUncontended" ${showUncontendedContigs ? 'checked' : ''}> Show uncontended contigs
      </label></div>
      <div class="row"><label title="Off by default: a connected MAG shows as a bare hub for its contended contigs only. Turn on to also show its entire contig set.">
        <input type="checkbox" id="toggleConnectedMags" ${showConnectedMagContigs ? 'checked' : ''}> Show contigs for connected MAGs
      </label></div>
      <div class="row"><label>Arrange</label><select id="networkAlgorithm">
        <option value="ring">Ring (selected MAG centred)</option>
        <option value="petal">Petal (grouped by MAG)</option>
        <option value="force">Force-directed</option>
      </select></div>
      <div class="row-count" id="reconciliationNetworkNote"></div>
      <div id="reconciliationNetwork"></div>
      <div id="contigEvidence"></div>
    </div>
  `;
}

/**
 * The primitive the redesigned network is built around: given a selected
 * MAG, every *contended* contig any tool voted for it — not just the ones
 * bin-reconciliation.js's majority-vote bookkeeping happened to assign it
 * (mag.coreContigIds/disputedContigIds), which by construction excludes a
 * genuinely tied contig (see bin-reconciliation.js's majorityMagId
 * handling) from every MAG's list. Using the vote-based reverse index
 * (latest.contigIdsByMagId) instead means a tied contig still shows up
 * here, attached to every MAG it's tied between, rather than being
 * invisible from all of them. Plus, for every *other* MAG any tool voted
 * for instead of the selected one, on any of those contended contigs, that
 * MAG appears as a hub too, so the contested vote is visible from both
 * sides. One hop only, off the selected MAG's own contigs; a secondary
 * MAG's unrelated disputes/ties with a third MAG are not pulled in, so the
 * neighborhood stays a bounded, readable "this MAG and what it's
 * contending with," not a recursive expansion across the whole
 * reconciliation graph.
 *
 * Two things are deliberately left out unless asked for, via `opts`,
 * since pulling them in by default buried the actual dispute under
 * everything else nearby:
 * - `showUncontendedContigs` (default false): the selected MAG's own
 *   contigs that no tool disagreed on at all. Nothing to resolve there.
 * - `showConnectedMagContigs` (default false): a secondary MAG's *entire*
 *   contig set, not just the contended ones it shares with the selected
 *   MAG. Off by default, a secondary MAG shows as a bare hub — present so
 *   the contended vote has somewhere to point, without dragging in every
 *   other contig that MAG happens to contain.
 *
 * A vote can name a MAG that lost every contig's majority vote and so was
 * dropped from result.putativeMags entirely (bin-reconciliation.js's
 * nonEmptyMags filter, documented there) — filtered out via magsById.has
 * everywhere below, so such a MAG never becomes a hub with no contigs to
 * back it, and edges never point at a hub that doesn't exist.
 *
 * @param {{showUncontendedContigs?:boolean, showConnectedMagContigs?:boolean}} [opts]
 * @returns {{hubs, leaves, edges, truncated, totalContigs, allTotalContigs,
 *   secondaryCount, hiddenUncontendedCount, hiddenConnectedCount}|null}
 *   null if selectedMagId no longer names a real MAG
 */
function buildMagNeighborhood(selectedMagId, result, opts = {}) {
  const includeUncontended = opts.showUncontendedContigs || false;
  const includeConnected = opts.showConnectedMagContigs || false;

  const magsById = new Map(result.putativeMags.map((m) => [m.magId, m]));
  if (!magsById.has(selectedMagId)) return null;

  const primaryContigIds = latest.contigIdsByMagId.get(selectedMagId) || new Set();

  // Weight each secondary MAG by how many of the selected MAG's contigs
  // actually contend with it. Real multi-tool data can have a MAG
  // "contending" with hundreds of secondaries — a fragmenting tool (VAMB
  // producing near-one-bin-per-contig output, on a real dataset that
  // surfaced this) means the overwhelming majority of those secondaries
  // share exactly one contig with the selected MAG, pure noise next to
  // the handful that are a genuine, substantial rival. Showing every one
  // as its own hub, unbounded, makes the ring unreadable and buries the
  // real dispute — only the top HUB_LIMIT by weight become hubs.
  const secondaryWeight = new Map();
  for (const contigId of primaryContigIds) {
    const entry = latest.contigAgreementEntryByContigId.get(contigId);
    if (!entry || entry.distinctGroupsVoted <= 1) continue;
    for (const magId of Object.values(entry.votes)) {
      if (magId && magId !== selectedMagId && magsById.has(magId)) {
        secondaryWeight.set(magId, (secondaryWeight.get(magId) || 0) + 1);
      }
    }
  }
  const HUB_LIMIT = 40;
  const allSecondaryMagIds = [...secondaryWeight.keys()];
  const rankedSecondaryMagIds = allSecondaryMagIds.sort((a, b) => secondaryWeight.get(b) - secondaryWeight.get(a) || a.localeCompare(b));
  const secondaryMagIds = new Set(rankedSecondaryMagIds.slice(0, HUB_LIMIT));
  const hiddenSecondaryCount = allSecondaryMagIds.length - secondaryMagIds.size;

  const shownMagIds = new Set([selectedMagId, ...secondaryMagIds]);
  const shownMags = [...shownMagIds].map((id) => magsById.get(id));

  // The selected MAG's own contigs: contended ones (entry.distinctGroupsVoted
  // > 1 — bin-reconciliation.js's count of how many distinct MAGs got a vote
  // for this contig) show only if at least one of the *other* MAGs they
  // contend with actually made the hub cut above — a contig whose sole
  // rival was trimmed as noise is left out of this view entirely (counted
  // in hiddenLowRelevanceCount) rather than drawn as if it were
  // uncontested, which would misrepresent it. Uncontended contigs only
  // show when that toggle is on.
  const contigIdSet = new Set();
  let hiddenUncontendedCount = 0;
  let hiddenLowRelevanceCount = 0;
  for (const contigId of primaryContigIds) {
    const entry = latest.contigAgreementEntryByContigId.get(contigId);
    const contended = entry && entry.distinctGroupsVoted > 1;
    if (!contended) {
      if (includeUncontended) contigIdSet.add(contigId);
      else hiddenUncontendedCount++;
      continue;
    }
    const stillContendedInView = Object.values(entry.votes).some((magId) => magId && magId !== selectedMagId && shownMagIds.has(magId));
    if (stillContendedInView) contigIdSet.add(contigId);
    else hiddenLowRelevanceCount++;
  }

  // A secondary MAG's own full contig set only shows when that toggle is
  // on — its contended contigs are already included above via the
  // selected MAG's own scan, so this only ever adds contigs unrelated to
  // the actual dispute.
  let hiddenConnectedCount = 0;
  for (const magId of secondaryMagIds) {
    const magContigIds = latest.contigIdsByMagId.get(magId) || new Set();
    for (const id of magContigIds) {
      if (contigIdSet.has(id)) continue;
      if (includeConnected) contigIdSet.add(id);
      else hiddenConnectedCount++;
    }
  }

  const allContigIds = [...contigIdSet];

  const NODE_LIMIT = 300;
  const contigIds = allContigIds.slice(0, NODE_LIMIT);
  const truncated = allContigIds.length > NODE_LIMIT;

  // A leaf's state drives the network's colouring (reconciliation-network.js):
  // 'core' (only ever voted into one shown MAG, nothing to decide), 'tied'
  // (2+ shown MAGs, no single majority winner among the tools that voted —
  // bin-reconciliation.js's majorityMagId is null, so there is no default
  // to fall back on at all), 'disputed' (2+ shown MAGs, one has a clear
  // majority but the student hasn't confirmed it), 'resolved' (2+ shown
  // MAGs, the student has explicitly assigned it via the evidence panel),
  // 'excluded' (the student explicitly removed it from both/all).
  // "Explicitly" matters for resolved/excluded: workingAssignment already
  // carries a majority-vote default for every non-tied voted contig from
  // deriveInitialAssignment, so decidedContigIds — populated only by the
  // evidence panel's own buttons — is what actually distinguishes "still
  // needs a look" from "you decided this," not merely "workingAssignment
  // has an entry."
  const leaves = [];
  const edges = [];
  for (const contigId of contigIds) {
    const entry = latest.contigAgreementEntryByContigId.get(contigId);
    if (!entry) continue;
    const hubIds = [...new Set(Object.values(entry.votes).filter((id) => id && shownMagIds.has(id)))];
    if (hubIds.length === 0) continue;
    const decided = decidedContigIds.has(contigId);
    let state;
    if (decided && workingAssignment.get(contigId) === EXCLUDED_BIN_ID) state = 'excluded';
    else if (hubIds.length === 1) state = 'core';
    else if (decided) state = 'resolved';
    else if (entry.majorityMagId === null) state = 'tied';
    else state = 'disputed';
    leaves.push({ id: contigId, hubIds, state });
    for (const [tool, magId] of Object.entries(entry.votes)) {
      if (magId && shownMagIds.has(magId)) edges.push({ leafId: contigId, hubId: magId, tool });
    }
  }

  const hubs = shownMags.map((m) => ({ id: m.magId, label: m.magId }));
  return {
    hubs, leaves, edges, truncated,
    totalContigs: contigIds.length, allTotalContigs: allContigIds.length,
    secondaryCount: secondaryMagIds.size, totalSecondaryCount: allSecondaryMagIds.length, hiddenSecondaryCount,
    tiedCount: leaves.filter((l) => l.state === 'tied').length,
    hiddenUncontendedCount, hiddenConnectedCount, hiddenLowRelevanceCount,
  };
}

/**
 * Renders the hub-and-leaf network into the placeholder left by
 * renderReconciliationCard, above — split out because it needs a live DOM
 * node (createReconciliationNetwork builds SVG into it), unlike the rest
 * of that card which is a plain innerHTML string. `neighborhood` is
 * `null` until a MAG is selected. The chosen layout algorithm
 * (`networkAlgorithm`, module-level) persists across re-renders and
 * changing it triggers a full re-render (simplest way to keep the network
 * and the rest of the page — the evidence panel especially — consistent).
 */
function initMagNetwork(result, neighborhood, records) {
  const { createReconciliationNetwork } = window.ClannMAG.reconciliationNetwork;
  const container = document.getElementById('reconciliationNetwork');
  const note = document.getElementById('reconciliationNetworkNote');
  const evidenceContainer = document.getElementById('contigEvidence');
  const algorithmSelect = document.getElementById('networkAlgorithm');
  const uncontendedToggle = document.getElementById('toggleUncontended');
  const connectedToggle = document.getElementById('toggleConnectedMags');
  if (!container || !note) return;

  // These three controls are rebuilt fresh (innerHTML) on every render, so
  // there's nothing stale to guard against re-wiring — unlike a control
  // that persists across renders, each one only ever gets wired once per
  // element instance. All three trigger a full re-render on change, same
  // as every other state toggle in this app (filters, params, MAG
  // selection), so the network, the note text, and the evidence panel all
  // stay consistent with whatever scope is now in effect.
  if (algorithmSelect) {
    algorithmSelect.value = networkAlgorithm;
    algorithmSelect.addEventListener('change', () => {
      networkAlgorithm = algorithmSelect.value;
      renderFilteredExplorer();
    });
  }
  if (uncontendedToggle) {
    uncontendedToggle.addEventListener('change', () => {
      showUncontendedContigs = uncontendedToggle.checked;
      renderFilteredExplorer();
    });
  }
  if (connectedToggle) {
    connectedToggle.addEventListener('change', () => {
      showConnectedMagContigs = connectedToggle.checked;
      renderFilteredExplorer();
    });
  }

  if (!selectedMagId || !neighborhood) {
    container.innerHTML = '';
    note.textContent = 'Select a putative MAG above to see its contended contigs and any MAGs it contends with.';
    if (evidenceContainer) evidenceContainer.innerHTML = '';
    return;
  }

  note.textContent = `${selectedMagId}: ${neighborhood.totalContigs.toLocaleString()} contig(s) shown` +
    (neighborhood.totalSecondaryCount > 0
      ? ` across ${neighborhood.hubs.length} MAG(s) (contending with ${neighborhood.totalSecondaryCount.toLocaleString()})`
      : ', no contested overlaps with other MAGs') +
    (neighborhood.hiddenSecondaryCount > 0 ? ` — ${neighborhood.hiddenSecondaryCount.toLocaleString()} more rival MAG(s) share only a contig or two each, hidden for readability` : '') +
    (neighborhood.tiedCount > 0 ? ` — ${neighborhood.tiedCount} tied, with no majority vote at all (bright ring)` : '') +
    (neighborhood.hiddenUncontendedCount > 0 ? ` — ${neighborhood.hiddenUncontendedCount.toLocaleString()} uncontended contig(s) hidden, toggle above to show` : '') +
    (neighborhood.hiddenConnectedCount > 0 ? ` — ${neighborhood.hiddenConnectedCount.toLocaleString()} more contig(s) available from connected MAGs, toggle above to show` : '') +
    (neighborhood.hiddenLowRelevanceCount > 0 ? ` — ${neighborhood.hiddenLowRelevanceCount.toLocaleString()} disputed contig(s) hidden (their only rival MAG isn't shown)` : '') +
    (neighborhood.truncated ? ` — truncated from ${neighborhood.allTotalContigs.toLocaleString()} for readability, narrow with MAG filters` : '') +
    '. Click a contig to compare evidence and decide where it belongs (selection stays highlighted); click another MAG to explore its neighborhood.';

  createReconciliationNetwork(container, { hubs: neighborhood.hubs, leaves: neighborhood.leaves, edges: neighborhood.edges }, {
    width: 680, height: 680, algorithm: networkAlgorithm, centralHubId: selectedMagId, selectedLeafId: selectedContigId,
    onLeafClick: (leafId) => { selectedContigId = leafId; renderContigEvidence(leafId, result, records); },
    onHubClick: (hubId) => { selectedMagId = hubId; selectedContigId = null; renderFilteredExplorer(); },
  });

  if (selectedContigId) renderContigEvidence(selectedContigId, result, records);
}

/**
 * The click-to-decide evidence panel (this simplification's central new
 * feature): for a clicked contig, one row per MAG any tool voted it into,
 * each with the evidence a student would weigh by hand — how this
 * contig's GC%/coverage compares to that MAG's other core contigs, and
 * whether each of its marker-gene families would be a unique contribution
 * there or a duplicate of one already present (computeMarkerContributions,
 * reused from the outlier-flagging module, run hypothetically against
 * "this MAG's core set plus this one contig" rather than the contig's
 * actual current bin). "Assign here" and "Exclude from both/all" both
 * just call working-assignment.js's one reassignment primitive — the
 * network, the MAG picker table, and this panel all read the same
 * `workingAssignment` Map, so a decision here is immediately visible in
 * all three on the next render.
 */
function renderContigEvidence(contigId, result, records) {
  const container = document.getElementById('contigEvidence');
  if (!container) return;

  const recordsById = new Map(records.map((r) => [r.id, r]));
  const record = recordsById.get(contigId);
  const entry = latest.contigAgreementEntryByContigId.get(contigId);
  if (!record || !entry) {
    container.innerHTML = `<div class="hint">No evidence available for ${contigId}.</div>`;
    return;
  }

  const { computeMarkerContributions } = window.ClannMAG.binSummary;
  const { reassignContigs } = window.ClannMAG.workingAssignment;
  const magsById = new Map(result.putativeMags.map((m) => [m.magId, m]));
  const tools = result.tools;

  const candidateMagIds = [...new Set(Object.values(entry.votes).filter(Boolean))].filter((id) => magsById.has(id));
  const currentDecision = workingAssignment.get(contigId) || null;

  const contigGc = record.gcContent * 100;
  const contigCov = record.coverageDepths ? record.coverageDepths.reduce((a, b) => a + b, 0) / record.coverageDepths.length : null;
  const contigFamilies = [...new Set((record.markerHits || []).map((h) => h.family))];

  const rows = candidateMagIds
    .map((magId) => {
      const mag = magsById.get(magId);
      const coreRecords = mag.coreContigIds
        .map((id) => recordsById.get(id))
        .filter((r) => r && r.id !== contigId);

      const meanGc = coreRecords.length ? (coreRecords.reduce((s, r) => s + r.gcContent, 0) / coreRecords.length) * 100 : null;
      const covRecords = coreRecords.filter((r) => r.coverageDepths);
      const meanCov = covRecords.length
        ? covRecords.reduce((s, r) => s + r.coverageDepths.reduce((a, b) => a + b, 0) / r.coverageDepths.length, 0) / covRecords.length
        : null;
      const gcDiff = meanGc !== null ? contigGc - meanGc : null;
      const covRatio = meanCov !== null && contigCov !== null && meanCov > 0 ? contigCov / meanCov : null;

      // "If this contig were added to this MAG's core set, would each of
      // its marker families be a new (unique) one, or a duplicate
      // (redundant) of a family already present elsewhere in the MAG?" —
      // answered by running the same unique/redundant logic the outlier
      // card uses, on a synthetic "core set + this contig" grouping.
      const contributions = computeMarkerContributions([...coreRecords, record]).get(contigId);
      const markerStatus = contigFamilies.length === 0
        ? '<span class="hint">none</span>'
        : contigFamilies
          .map((fam) => {
            const isRedundant = contributions.redundantFamilies.includes(fam);
            return `<span class="marker-tag ${isRedundant ? 'marker-redundant' : 'marker-unique'}" title="${isRedundant ? 'already present elsewhere in this MAG — assigning here would duplicate it' : 'not currently found elsewhere in this MAG — assigning here would raise completeness'}">${fam}</span>`;
          })
          .join(' ');

      // One column per loaded tool showing *that tool's own bin ID* if its
      // vote for this contig landed in this MAG, blank otherwise — not a
      // flattened text list (mag.members already carries the tool->binId
      // mapping; entry.votes says which MAG each tool actually voted for,
      // so a blank cell here means that tool voted elsewhere or had no
      // opinion, not that data is missing).
      const toolCells = tools
        .map((tool) => {
          if (entry.votes[tool] !== magId) return '<td class="hint">—</td>';
          const member = mag.members.find((m) => m.tool === tool);
          return `<td>${member ? member.binId : ''}</td>`;
        })
        .join('');

      const isCurrent = currentDecision === magId;
      return `<tr class="${isCurrent ? 'evidence-row-current' : ''}">
        <td>${magId}</td>
        ${toolCells}
        <td class="num">${meanGc === null ? '<span class="hint">n/a</span>' : `${meanGc.toFixed(1)}%`}</td>
        <td class="num">${gcDiff === null ? '<span class="hint">n/a</span>' : `${gcDiff > 0 ? '+' : ''}${gcDiff.toFixed(1)}pp`}</td>
        <td class="num">${meanCov === null ? '<span class="hint">n/a</span>' : meanCov.toFixed(1)}</td>
        <td class="num">${covRatio === null ? '<span class="hint">n/a</span>' : `${covRatio.toFixed(2)}&times;`}</td>
        <td>${markerStatus}</td>
        <td>${isCurrent ? '<strong>current</strong>' : `<button class="act evidence-assign-btn" type="button" data-mag-id="${magId}">Assign here</button>`}</td>
      </tr>`;
    })
    .join('');

  const disputeNote = candidateMagIds.length <= 1
    ? 'not contested — every voting tool agrees'
    : entry.majorityMagId === null
      ? `tied between ${candidateMagIds.length} MAGs across ${entry.totalVotes} tool vote(s) — no majority, this one's your call`
      : `contested between ${candidateMagIds.length} MAGs across ${entry.totalVotes} tool vote(s) (${entry.majorityMagId} currently leads the vote)`;

  container.innerHTML = `
    <div class="card evidence-card">
      <h3>Contig ${contigId}</h3>
      <div class="row-count">Length ${record.length.toLocaleString()} bp &middot; own GC ${contigGc.toFixed(1)}%${contigCov !== null ? ` &middot; own mean coverage ${contigCov.toFixed(1)}` : ''} &middot; ${disputeNote}</div>
      <div class="table-wrap scroll-panel">
        <table class="data-table">
          <thead><tr>
            <th>Candidate MAG</th>
            ${tools.map((t) => `<th title="This tool's own bin ID, if it voted this contig into this MAG">${t}</th>`).join('')}
            <th title="Mean GC% of this MAG's other core (unanimous-agreement) contigs">MAG mean GC</th>
            <th title="This contig's GC% minus the MAG's mean core GC%">GC diff</th>
            <th title="Mean coverage depth of this MAG's other core contigs">MAG mean cov.</th>
            <th title="This contig's mean coverage divided by the MAG's mean core coverage">Cov. ratio</th>
            <th title="This contig's marker-gene families, and whether assigning it here would be a unique contribution or a duplicate">Marker genes</th>
            <th>Decision</th>
          </tr></thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
      <div class="row" style="margin-top:8px">
        <button class="act warn" id="evidenceExcludeBtn" type="button">${currentDecision === EXCLUDED_BIN_ID ? 'Excluded ✓' : 'Exclude from both/all'}</button>
      </div>
    </div>
  `;

  container.querySelectorAll('.evidence-assign-btn').forEach((btn) => {
    btn.addEventListener('click', () => {
      workingAssignment = reassignContigs(workingAssignment, [contigId], btn.dataset.magId);
      decidedContigIds.add(contigId);
      renderFilteredExplorer();
    });
  });
  const excludeBtn = document.getElementById('evidenceExcludeBtn');
  if (excludeBtn) {
    excludeBtn.addEventListener('click', () => {
      workingAssignment = reassignContigs(workingAssignment, [contigId], EXCLUDED_BIN_ID);
      decidedContigIds.add(contigId);
      renderFilteredExplorer();
    });
  }
}

/**
 * The working assignment (working-assignment.js) is the one live-editable
 * session state everything in this redesign reads/writes: the MAG picker
 * table's live numbers, the network's resolved/excluded leaf colouring,
 * the evidence panel's decisions, and Export. Initialized once per load
 * (guarded by workingAssignmentInitialized, same as before this
 * redesign — this used to run lazily inside the old Phase 7 scatter
 * section, moved here since that section is gone but the same
 * once-per-load initialization still needs to happen before anything
 * reads workingAssignment).
 */
function ensureWorkingAssignment() {
  if (workingAssignmentInitialized) return;
  const { deriveInitialAssignment } = window.ClannMAG.workingAssignment;
  workingAssignment = deriveInitialAssignment(latest.binTablesByTool, latest.reconciliationResult);
  workingAssignmentInitialized = true;
}

/** Triggers a browser save of `blob` as `filename` via a throwaway object URL. */
function triggerDownload(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/**
 * Slices a MAG/bin's contigs' raw sequence bytes straight out of the
 * original assembly File (brief §Export — "built via Blob.slice against
 * the first-pass index rather than a second full file read") and
 * reassembles them under their original headers into one multi-FASTA
 * Blob. A gzip source needs one full decompression pass first (Phase 1's
 * documented limitation: `.fai`-style offsets are in the *decompressed*
 * stream) — done at most once per export, not per contig.
 * @returns {Promise<{blob: Blob, skippedContigIds: string[]}>}
 */
async function extractBinFasta(binId, contigIds) {
  const { planFastaExtraction } = window.ClannMAG.fastaExtract;
  const plan = planFastaExtraction(currentRecords, new Map([[binId, contigIds]]));
  const { entries, skippedContigIds } = plan.get(binId);

  const needsDecompression = entries.some((e) => {
    const record = currentRecords.find((r) => r.id === e.id);
    return record && record.faiEntry.sourceCompressed;
  });
  let decompressedBytes = null;
  if (needsDecompression) {
    const stream = currentAssemblyFile.stream().pipeThrough(new DecompressionStream('gzip'));
    decompressedBytes = new Uint8Array(await new Response(stream).arrayBuffer());
  }

  const parts = [];
  for (const entry of entries) {
    parts.push(`>${entry.header}\n`);
    const { offset, byteLength } = entry.span;
    const seqBytes = decompressedBytes
      ? decompressedBytes.slice(offset, offset + byteLength)
      : new Uint8Array(await currentAssemblyFile.slice(offset, offset + byteLength).arrayBuffer());
    parts.push(seqBytes);
    if (seqBytes.length === 0 || seqBytes[seqBytes.length - 1] !== 10) parts.push('\n');
  }
  return { blob: new Blob(parts, { type: 'text/plain' }), skippedContigIds };
}

/**
 * Phase 9: export section — revised contig->bin assignment table (CSV),
 * per-MAG/bin summary table (CSV), and per-bin FASTA extraction. Reads
 * from the current working assignment, so an export always reflects any
 * manual reassignment made in the session (brief's "revised" framing),
 * not the originally loaded tables.
 */
function initExportSection() {
  const { assignmentToRows, listBinIds } = window.ClannMAG.workingAssignment;
  const { computeBinSummaries } = window.ClannMAG.binSummary;
  const { assignmentToCsv, binSummaryToCsv } = window.ClannMAG.exportCsv;

  const card = document.getElementById('export-card');
  card.innerHTML = `
    <h3>Export</h3>
    <div class="row-count">Exports reflect the current working assignment, including any reassignments made above.</div>
    <div class="row"><button class="act" id="exportAssignmentBtn" type="button">Download revised assignment table (CSV)</button></div>
    <div class="row"><button class="act" id="exportSummaryBtn" type="button">Download bin summary table (CSV)</button></div>
    <div class="row">
      <select id="exportBinSelect"></select>
      <button class="act" id="exportFastaBtn" type="button">Download bin FASTA</button>
    </div>
    <div class="row-count" id="exportNote"></div>
  `;

  function refreshBinOptions() {
    document.getElementById('exportBinSelect').innerHTML =
      listBinIds(workingAssignment).map((b) => `<option value="${b}">${b}</option>`).join('');
  }
  refreshBinOptions();
  document.getElementById('exportBinSelect').addEventListener('focus', refreshBinOptions);

  document.getElementById('exportAssignmentBtn').addEventListener('click', () => {
    const csv = assignmentToCsv(assignmentToRows(workingAssignment));
    triggerDownload(new Blob([csv], { type: 'text/csv' }), 'contig-bin-assignment.csv');
  });

  document.getElementById('exportSummaryBtn').addEventListener('click', () => {
    refreshBinOptions();
    const { summaries } = computeBinSummaries(currentRecords, assignmentToRows(workingAssignment), { thresholds: currentParams.mimag, recallRate: currentParams.recallRate });
    const csv = binSummaryToCsv(summaries);
    triggerDownload(new Blob([csv], { type: 'text/csv' }), 'bin-summary.csv');
  });

  document.getElementById('exportFastaBtn').addEventListener('click', async () => {
    const binId = document.getElementById('exportBinSelect').value;
    if (!binId) return;
    const note = document.getElementById('exportNote');
    note.textContent = `Extracting ${binId}…`;
    const contigIds = [...workingAssignment.entries()].filter(([, b]) => b === binId).map(([c]) => c);
    const { blob, skippedContigIds } = await extractBinFasta(binId, contigIds);
    triggerDownload(blob, `${binId}.fasta`);
    note.textContent = skippedContigIds.length > 0
      ? `${skippedContigIds.length.toLocaleString()} contig(s) skipped (non-uniform FASTA line wrapping can't be safely re-sliced).`
      : '';
  });
}

/**
 * Renders the left-pane filters UI (brief §Left pane: filters and search)
 * into #filtersSectionBody, against whatever's in `latest`. Rebuilt
 * wholesale on every fresh load (bin options depend on which tools were
 * loaded), not on every filter change — inputs keep their own DOM state
 * between keystrokes, `currentFilters` is only read from them on change.
 */
function renderFiltersSection() {
  const { listBinFilterOptions } = window.ClannMAG.filters;
  const body = document.getElementById('filtersSectionBody');
  const binOptions = listBinFilterOptions(latest.binIndex);

  body.innerHTML = `
    <div class="row"><label>Search</label><input type="text" id="filterSearch" placeholder="contig or bin ID"></div>
    <div class="row"><label>Length ≥ bp</label><input type="number" id="filterLengthMin" min="0"></div>
    <div class="row"><label>Length ≤ bp</label><input type="number" id="filterLengthMax" min="0"></div>
    <div class="row"><label>GC% ≥</label><input type="number" id="filterGcMin" min="0" max="100"></div>
    <div class="row"><label>GC% ≤</label><input type="number" id="filterGcMax" min="0" max="100"></div>
    <div class="row"><label>Coding density% ≥</label><input type="number" id="filterCdMin" min="0" max="100"></div>
    <div class="row"><label>Coding density% ≤</label><input type="number" id="filterCdMax" min="0" max="100"></div>
    ${binOptions.length > 0 ? `
    <div class="row"><label>Bin</label><select id="filterBin">
      <option value="">Any</option>
      <option value="__unbinned__">Unbinned (all tools)</option>
      ${binOptions.map((o) => `<option value="${o.value}">${o.label}</option>`).join('')}
    </select></div>` : ''}
    ${latest.reconciliationResult ? `
    <div class="row"><label title="Show only contigs where the loaded tools agree on the assigned MAG in at most this fraction of votes">Max agreement %</label>
      <input type="number" id="filterMaxAgreement" min="0" max="100" placeholder="100"></div>` : ''}
    <div class="row"><button class="act" id="filterResetBtn" type="button">Reset filters</button></div>
    <div class="hint" id="filterSummary"></div>
  `;

  function readFiltersFromInputs() {
    const num = (id) => {
      const el = document.getElementById(id);
      if (!el || el.value === '') return null;
      const v = Number(el.value);
      return Number.isNaN(v) ? null : v;
    };
    currentFilters = {
      lengthMin: num('filterLengthMin'), lengthMax: num('filterLengthMax'),
      gcMin: num('filterGcMin'), gcMax: num('filterGcMax'),
      codingDensityMin: num('filterCdMin'), codingDensityMax: num('filterCdMax'),
      binFilter: document.getElementById('filterBin')?.value || '',
      maxAgreementPercent: num('filterMaxAgreement'),
      searchText: document.getElementById('filterSearch')?.value || '',
    };
    renderFilteredExplorer();
  }

  // Debounced so typing in the search box (or a number field) doesn't
  // re-render every keystroke — the debounce is short enough to still
  // feel live.
  let debounceHandle = null;
  const onInput = () => {
    clearTimeout(debounceHandle);
    debounceHandle = setTimeout(readFiltersFromInputs, 150);
  };
  body.querySelectorAll('input').forEach((el) => el.addEventListener('input', onInput));
  body.querySelectorAll('select').forEach((el) => el.addEventListener('change', readFiltersFromInputs));
  document.getElementById('filterResetBtn').addEventListener('click', () => {
    currentFilters = window.ClannMAG.filters.defaultFilters();
    renderFiltersSection();
    renderFilteredExplorer();
  });
}

/**
 * Renders the MAG-filters UI (a second, independent filter axis over the
 * cross-tool reconciliation table's own columns — see currentMagFilters
 * above) into #magFiltersSectionBody. Only meaningful once 2+ bin tables
 * are loaded and reconciliation has run; otherwise shows a hint instead.
 * Rebuilt on every fresh load (the "supported by tool" dropdown depends
 * on which tools were loaded), not on every filter change.
 */
function renderMagFiltersSection() {
  const body = document.getElementById('magFiltersSectionBody');
  if (!latest.reconciliationResult) {
    body.innerHTML = '<div class="hint">Load 2+ contig&rarr;bin tables to enable cross-tool reconciliation, then filter by MAG here.</div>';
    return;
  }
  const tools = latest.tools;

  body.innerHTML = `
    <div class="row"><label>MAG ID</label><input type="text" id="magFilterSearch" placeholder="e.g. MAG_1"></div>
    <div class="row"><label>Tier</label></div>
    <div class="row"><label><input type="checkbox" id="magFilterTierHigh" checked> High</label></div>
    <div class="row"><label><input type="checkbox" id="magFilterTierMedium" checked> Medium</label></div>
    <div class="row"><label><input type="checkbox" id="magFilterTierLow" checked> Low</label></div>
    <div class="row"><label>Core ≥</label><input type="number" id="magFilterCoreMin" min="0"></div>
    <div class="row"><label>Core ≤</label><input type="number" id="magFilterCoreMax" min="0"></div>
    <div class="row"><label>Disputed ≥</label><input type="number" id="magFilterDisputedMin" min="0"></div>
    <div class="row"><label>Disputed ≤</label><input type="number" id="magFilterDisputedMax" min="0"></div>
    <div class="row"><label>Completeness% ≥</label><input type="number" id="magFilterCompletenessMin" min="0" max="100"></div>
    <div class="row"><label>Completeness% ≤</label><input type="number" id="magFilterCompletenessMax" min="0" max="100"></div>
    <div class="row"><label>Redundancy% ≥</label><input type="number" id="magFilterRedundancyMin" min="0" max="100"></div>
    <div class="row"><label>Redundancy% ≤</label><input type="number" id="magFilterRedundancyMax" min="0" max="100"></div>
    <div class="row"><label>Supported by</label><select id="magFilterTool">
      <option value="">Any tool</option>
      ${tools.map((t) => `<option value="${t}">${t}</option>`).join('')}
    </select></div>
    <div class="row"><button class="act" id="magFilterResetBtn" type="button">Reset MAG filters</button></div>
    <div class="hint" id="magFilterSummary"></div>
  `;

  function readMagFiltersFromInputs() {
    const num = (id) => {
      const el = document.getElementById(id);
      if (!el || el.value === '') return null;
      const v = Number(el.value);
      return Number.isNaN(v) ? null : v;
    };
    currentMagFilters = {
      magIdSearch: document.getElementById('magFilterSearch')?.value || '',
      tiers: {
        high: document.getElementById('magFilterTierHigh').checked,
        medium: document.getElementById('magFilterTierMedium').checked,
        low: document.getElementById('magFilterTierLow').checked,
      },
      coreMin: num('magFilterCoreMin'), coreMax: num('magFilterCoreMax'),
      disputedMin: num('magFilterDisputedMin'), disputedMax: num('magFilterDisputedMax'),
      completenessMin: num('magFilterCompletenessMin'), completenessMax: num('magFilterCompletenessMax'),
      redundancyMin: num('magFilterRedundancyMin'), redundancyMax: num('magFilterRedundancyMax'),
      supportedByTool: document.getElementById('magFilterTool')?.value || '',
    };
    renderFilteredExplorer();
  }

  let debounceHandle = null;
  const onInput = () => {
    clearTimeout(debounceHandle);
    debounceHandle = setTimeout(readMagFiltersFromInputs, 150);
  };
  body.querySelectorAll('input[type=text], input[type=number]').forEach((el) => el.addEventListener('input', onInput));
  body.querySelectorAll('input[type=checkbox], select').forEach((el) => el.addEventListener('change', readMagFiltersFromInputs));
  document.getElementById('magFilterResetBtn').addEventListener('click', () => {
    currentMagFilters = window.ClannMAG.magFilters.defaultMagFilters();
    renderMagFiltersSection();
    renderFilteredExplorer();
  });
}

/**
 * Excludes/re-includes one tool from cross-tool reconciliation — the
 * mechanism the "Binning tools" section exposes for a tool whose output
 * looks unreliable (a real dataset surfaced this: a tool whose bins were
 * ~99% singletons dragged completeness/agreement toward zero everywhere
 * it had an opinion, purely by disagreeing with everyone else almost at
 * random). Refuses to drop the last active tool — reconciliation needs
 * at least one — rather than leaving every bin-derived view empty with
 * no obvious way back. Session state that depends on *which* MAGs exist
 * (selection, decisions, the working assignment) is reset, same as a
 * fresh load: which putative MAGs even exist can change completely once
 * a tool's votes are added or removed, so anything referencing the old
 * MAG IDs would silently point at the wrong (or a nonexistent) thing.
 */
function setToolActive(tool, active) {
  if (!active && activeTools.size <= 1 && activeTools.has(tool)) return;
  if (active) activeTools.add(tool); else activeTools.delete(tool);

  selectedMagId = null;
  selectedContigId = null;
  decidedContigIds = new Set();
  workingAssignmentInitialized = false;
  networkAlgorithm = 'ring';

  recomputeLatest(latest.records, filterActiveBinTables()).then(() => {
    renderToolsSection();
    renderMagFiltersSection();
    renderFiltersSection();
    renderFilteredExplorer();
  });
}

/**
 * Renders the "Binning tools" left-pane section (index.html's previously
 * unwired #toolsSectionBody stub) — per-tool bin/contig counts and,
 * specifically, each tool's singleton-bin fraction, the single number
 * that would have flagged a real fragmented-tool problem immediately
 * rather than needing to be root-caused through the network view three
 * layers downstream. Always lists every *loaded* tool (allBinTablesByTool),
 * not just the active ones, so excluding a tool doesn't remove the only
 * way to re-include it.
 */
function renderToolsSection() {
  const body = document.getElementById('toolsSectionBody');
  if (!body) return;
  if (!allBinTablesByTool || allBinTablesByTool.size === 0) {
    body.innerHTML = '<div class="hint">Load one or more contig&rarr;bin tables to see per-tool stats here.</div>';
    return;
  }

  const { isUnbinnedLabel } = window.ClannMAG.binReconciliation;
  const rows = [...allBinTablesByTool.entries()]
    .map(([tool, assignments]) => {
      const binSizes = new Map();
      for (const { binId } of assignments) {
        if (isUnbinnedLabel(binId)) continue;
        binSizes.set(binId, (binSizes.get(binId) || 0) + 1);
      }
      const binCount = binSizes.size;
      const singletonCount = [...binSizes.values()].filter((n) => n === 1).length;
      const singletonPct = binCount ? (singletonCount / binCount) * 100 : 0;
      const isActive = activeTools.has(tool);
      const warn = singletonPct >= 50
        ? ' <span class="hint" title="More than half of this tool\'s bins contain only one contig — its votes may add more noise than signal to the reconciliation">&#9888;</span>'
        : '';
      return `<tr class="${isActive ? '' : 'tool-row-inactive'}">
        <td><label><input type="checkbox" class="tool-active-checkbox" data-tool="${tool}" ${isActive ? 'checked' : ''}> ${tool}</label></td>
        <td class="num">${assignments.length.toLocaleString()}</td>
        <td class="num">${binCount.toLocaleString()}</td>
        <td class="num">${singletonPct.toFixed(0)}%${warn}</td>
      </tr>`;
    })
    .join('');

  body.innerHTML = `
    <div class="hint">Exclude a tool to see how the reconciliation changes without it — useful when one tool's output looks unusually fragmented (high singleton-bin %). At least one tool must stay active.</div>
    <table class="data-table">
      <thead><tr>
        <th>Tool</th>
        <th title="Total contig rows this tool assigned to a bin">Contigs</th>
        <th title="Distinct bins this tool produced">Bins</th>
        <th title="Fraction of this tool's bins containing exactly one contig — high values suggest fragmented, low-signal output">Singleton bins</th>
      </tr></thead>
      <tbody>${rows}</tbody>
    </table>
  `;

  body.querySelectorAll('.tool-active-checkbox').forEach((cb) => {
    cb.addEventListener('change', () => {
      const tool = cb.dataset.tool;
      if (!cb.checked && activeTools.size <= 1 && activeTools.has(tool)) {
        cb.checked = true; // refuse — reconciliation needs at least one active tool
        return;
      }
      setToolActive(tool, cb.checked);
    });
  });
}

/**
 * Renders the Thresholds & parameters UI into #paramsSectionBody: every
 * global cutoff a student might want to see the effect of moving (MIMAG
 * quality-tier thresholds, the cross-tool bin-matching overlap threshold,
 * and outlier-flagging thresholds). Changing minJaccard re-runs bin matching itself
 * (recomputeLatest); every other field only changes what gets displayed
 * from already-computed data, so those just re-render.
 */
function renderParamsSection() {
  const body = document.getElementById('paramsSectionBody');
  const p = currentParams;

  body.innerHTML = `
    <div class="hint">These control how the tool's own derived calls (quality tiers, bin matching, outlier flags, duplicate detection) are drawn — adjust them to see how sensitive the results are to where the lines sit.</div>
    <div class="row"><label title="The built-in marker-gene search is a fast, approximate heuristic, not a profile-HMM search — an independent BLAST verification (docs/scg-blast-verification.md) measured it recovers only about this fraction of genuinely-present marker genes, even after this app's own threshold tuning. Completeness/redundancy below are divided by this rate rather than read raw off the family-hit count, which would otherwise understate every genuinely-complete genome by roughly this same gap.">Assumed marker-search recall</label><input type="number" id="paramRecallRate" min="0.01" max="1" step="0.01" value="${p.recallRate}"></div>
    <div class="row"><label title="Completeness/contamination proxy only, and now corrected for the assumed recall rate above — see the bin summary note">MIMAG High: completeness &gt;</label><input type="number" id="paramMimagHighComp" min="0" max="100" value="${p.mimag.highMinCompleteness}"></div>
    <div class="row"><label>MIMAG High: contamination &lt;</label><input type="number" id="paramMimagHighCont" min="0" max="100" value="${p.mimag.highMaxContamination}"></div>
    <div class="row"><label>MIMAG Medium: completeness ≥</label><input type="number" id="paramMimagMedComp" min="0" max="100" value="${p.mimag.mediumMinCompleteness}"></div>
    <div class="row"><label>MIMAG Medium: contamination &lt;</label><input type="number" id="paramMimagMedCont" min="0" max="100" value="${p.mimag.mediumMaxContamination}"></div>
    <div class="row"><label title="Reciprocal-best-hit overlap required before two tools' bins are matched as the same putative MAG">Min bin-match Jaccard</label><input type="number" id="paramMinJaccard" min="0" max="1" step="0.01" value="${p.minJaccard}"></div>
    <div class="row"><label title="Composition/coverage z-score above which a contig counts as an outlier flag">Outlier Z threshold</label><input type="number" id="paramZThreshold" min="0" step="0.1" value="${p.outlier.zThreshold}"></div>
    <div class="row"><label title="Marker-gene taxonomic distance above which a contig counts as an outlier flag">Tax. distance threshold</label><input type="number" id="paramTaxDistance" min="0" step="1" value="${p.outlier.taxDistanceThreshold}"></div>
    <div class="row"><button class="act" id="paramResetBtn" type="button">Reset to defaults</button></div>
  `;

  function readParamsFromInputs() {
    const num = (id, fallback) => {
      const el = document.getElementById(id);
      if (!el || el.value === '') return fallback;
      const v = Number(el.value);
      return Number.isNaN(v) ? fallback : v;
    };
    const previousMinJaccard = currentParams.minJaccard;
    currentParams = {
      mimag: {
        highMinCompleteness: num('paramMimagHighComp', p.mimag.highMinCompleteness),
        highMaxContamination: num('paramMimagHighCont', p.mimag.highMaxContamination),
        mediumMinCompleteness: num('paramMimagMedComp', p.mimag.mediumMinCompleteness),
        mediumMaxContamination: num('paramMimagMedCont', p.mimag.mediumMaxContamination),
      },
      minJaccard: num('paramMinJaccard', currentParams.minJaccard),
      outlier: {
        zThreshold: num('paramZThreshold', currentParams.outlier.zThreshold),
        taxDistanceThreshold: num('paramTaxDistance', currentParams.outlier.taxDistanceThreshold),
      },
      recallRate: Math.max(0.01, Math.min(1, num('paramRecallRate', currentParams.recallRate))),
    };
    if (currentParams.minJaccard !== previousMinJaccard) {
      recomputeLatest(latest.records, latest.binTablesByTool).then(renderFilteredExplorer);
    } else {
      renderFilteredExplorer();
    }
  }

  let debounceHandle = null;
  const onInput = () => {
    clearTimeout(debounceHandle);
    debounceHandle = setTimeout(readParamsFromInputs, 250);
  };
  body.querySelectorAll('input').forEach((el) => el.addEventListener('input', onInput));
  document.getElementById('paramResetBtn').addEventListener('click', () => {
    currentParams = defaultGlobalParams();
    renderParamsSection();
    recomputeLatest(latest.records, latest.binTablesByTool).then(renderFilteredExplorer);
  });
}

/**
 * Builds the right-pane analysis sections from `latest` (everything
 * filter-independent, computed once per load). Re-run on every filter,
 * param, MAG-selection, or working-assignment change; cheap enough (array
 * filter + string templating over records already in memory, no
 * re-parsing or re-searching) to do synchronously each time.
 *
 * Scope: `currentFilters` (left-pane contig filters) only ever narrows the
 * collapsed "All contigs" raw table at the bottom — every other view here
 * (the MAG picker, the contig network, the outlier card, Export) is keyed
 * off `selectedMagId`/the working assignment instead, not the contig
 * filters, since those views work at MAG/session-state granularity, not
 * per-contig-property granularity. `currentMagFilters` (left-pane MAG
 * filters) narrows which MAGs appear in the picker table — its job, per
 * the redesign, is helping find a MAG to select.
 */
function renderFilteredExplorer() {
  const { records, tools, reconciliationResult, outlierFlags, outlierMeta, binIndex, agreementByContigId } = latest;
  ensureWorkingAssignment();

  const { applyFilters } = window.ClannMAG.filters;
  const filteredRecords = applyFilters(records, currentFilters, { binIndex, agreementByContigId });

  const filterSummaryEl = document.getElementById('filterSummary');
  if (filterSummaryEl) {
    filterSummaryEl.textContent =
      `${filteredRecords.length.toLocaleString()} of ${records.length.toLocaleString()} contigs match the current filters (applies to the raw contig table below).`;
  }

  const sorted = [...filteredRecords].sort((a, b) => b.length - a.length);
  const totalLength = records.reduce((sum, r) => sum + r.length, 0);
  const n50 = computeN50([...records].sort((a, b) => b.length - a.length).map((r) => r.length), totalLength);
  const meanGc = records.length ? records.reduce((sum, r) => sum + r.gcContent, 0) / records.length : 0;
  const contigsWithMarkers = records.filter((r) => r.markerHits && r.markerHits.length > 0).length;
  const distinctFamiliesHit = new Set(records.flatMap((r) => (r.markerHits || []).map((h) => h.family))).size;

  const explorer = document.getElementById('explorer');
  const rows = sorted
    .map((r) => `<tr>
      <td>${r.id}</td>
      <td class="num">${r.length.toLocaleString()}</td>
      <td class="num">${(r.gcContent * 100).toFixed(1)}%</td>
      <td class="num">${r.gcSkew.toFixed(3)}</td>
      <td class="num">${(r.codingDensity * 100).toFixed(1)}%</td>
      <td>${formatMarkerHits(r.markerHits)}</td>
      <td class="num">${r.faiEntry.uniform ? '' : '⚠'}</td>
    </tr>`)
    .join('');

  let magSummaryData = [];
  let filteredMagIds = new Set();
  if (reconciliationResult) {
    magSummaryData = computeMagSummaryData(records, reconciliationResult);
    const { applyMagFilters } = window.ClannMAG.magFilters;
    filteredMagIds = new Set(applyMagFilters(magSummaryData, currentMagFilters).map((m) => m.magId));
  }
  const reconciliationCard = reconciliationResult
    ? renderReconciliationCard(records, reconciliationResult, magSummaryData, filteredMagIds)
    : '';

  const neighborhood = (reconciliationResult && selectedMagId)
    ? buildMagNeighborhood(selectedMagId, reconciliationResult, { showUncontendedContigs, showConnectedMagContigs })
    : null;
  const neighborhoodContigIds = neighborhood ? new Set(neighborhood.leaves.map((l) => l.id)) : null;
  const rankedOutlierFlags = applyOutlierThresholds(outlierFlags, currentParams.outlier);
  const scopedOutlierFlags = neighborhoodContigIds ? rankedOutlierFlags.filter((f) => neighborhoodContigIds.has(f.contigId)) : [];
  const outlierCard = tools.length > 0 ? renderOutlierCard(scopedOutlierFlags, outlierMeta, selectedMagId) : '';

  explorer.innerHTML = `
    <div class="card">
      <h3>Assembly summary</h3>
      <div class="row"><label>Contigs</label><strong>${records.length.toLocaleString()}</strong></div>
      <div class="row"><label>Total length</label><strong>${totalLength.toLocaleString()} bp</strong></div>
      <div class="row"><label>N50</label><strong>${n50.toLocaleString()} bp</strong></div>
      <div class="row"><label>Mean GC</label><strong>${(meanGc * 100).toFixed(1)}%</strong></div>
      <div class="row"><label>Contigs with marker genes</label><strong>${contigsWithMarkers.toLocaleString()}</strong></div>
      <div class="row"><label>Distinct marker families hit</label><strong>${distinctFamiliesHit} / 40</strong></div>
    </div>
    ${reconciliationCard}
    ${outlierCard}
    ${tools.length > 0 ? '<div class="card" id="export-card"></div>' : ''}
    <div class="card">
      <details id="allContigsDetails">
        <summary><h3 style="display:inline-block;margin:0">All contigs (raw table)</h3></summary>
        <div class="row-count" style="margin-top:8px">${sorted.length.toLocaleString()} of ${records.length.toLocaleString()} contigs match the current contig filters, longest first</div>
        <div class="table-wrap scroll-panel">
          <table class="data-table">
            <thead><tr><th>Contig</th><th>Length</th><th>GC%</th><th>GC skew</th><th>Coding density</th><th>Marker genes</th><th title="Non-uniform FASTA line wrapping">⚠</th></tr></thead>
            <tbody>${rows}</tbody>
          </table>
        </div>
      </details>
    </div>
  `;

  if (reconciliationResult) initMagNetwork(reconciliationResult, neighborhood, records);
  if (tools.length > 0) initExportSection();
}

/**
 * Builds everything filtering must NOT redo — cross-tool reconciliation
 * (using the current `minJaccard`), the combined outlier/disagreement
 * flags (marker-gene taxonomy fetch included, cached after the first
 * call), and the bin-index/agreement lookups filters.js needs — and
 * stores it all in `latest`. Called once on a fresh load, and again
 * whenever `minJaccard` changes in the Thresholds & parameters panel,
 * since that's the one adjustable parameter that changes bin *matching*
 * itself (which MAGs exist at all), not just how already-computed numbers
 * get displayed — see the other parameters' inline recompute in
 * renderFilteredExplorer for the cheaper case.
 */
async function recomputeLatest(records, binTablesByTool) {
  const tools = binTablesByTool ? [...binTablesByTool.keys()] : [];

  let reconciliationResult = null;
  if (tools.length > 1) {
    reconciliationResult = window.ClannMAG.binReconciliation.reconcileBins(binTablesByTool, { minJaccard: currentParams.minJaccard });
  }

  let outlierFlags = [];
  let outlierMeta = { hasCoverage: false, hasTaxonomy: false, hasKraken: false, hasCrossTool: false };
  if (tools.length > 0) {
    const tree = await getTaxonomyTree();
    outlierFlags = computeOutlierFlags(records, binTablesByTool, reconciliationResult, tree);
    outlierMeta = {
      hasCoverage: records.some((r) => r.coverageDepths),
      hasTaxonomy: tree !== null,
      hasKraken: records.some((r) => r.krakenTaxId != null),
      hasCrossTool: reconciliationResult !== null,
    };
  }

  const { buildBinIndex } = window.ClannMAG.filters;
  const binIndex = buildBinIndex(binTablesByTool);
  const agreementByContigId = new Map();
  const contigAgreementEntryByContigId = new Map();
  // Reverse index: every contig that got *any* tool's vote for a given MAG,
  // not just the ones bin-reconciliation.js's majority-vote bookkeeping
  // happened to assign it (mag.coreContigIds/disputedContigIds) — that
  // bookkeeping leaves a genuinely tied contig (no single majority winner)
  // out of every MAG's list entirely, so buildMagNeighborhood needs this
  // broader index to still find and display tied contigs when a MAG that
  // was one of the tied contenders gets selected.
  const contigIdsByMagId = new Map();
  if (reconciliationResult) {
    for (const c of reconciliationResult.contigAgreement) {
      agreementByContigId.set(c.contigId, c.agreementFraction);
      contigAgreementEntryByContigId.set(c.contigId, c);
      for (const magId of new Set(Object.values(c.votes).filter(Boolean))) {
        if (!contigIdsByMagId.has(magId)) contigIdsByMagId.set(magId, new Set());
        contigIdsByMagId.get(magId).add(c.contigId);
      }
    }
  }

  latest = {
    records, binTablesByTool, tools, reconciliationResult, outlierFlags, outlierMeta, binIndex,
    agreementByContigId, contigAgreementEntryByContigId, contigIdsByMagId,
  };
}

/**
 * Compute-once entry point for a freshly loaded assembly (+ any bin/
 * coverage/Kraken tables): resets every adjustable-parameter/filter state
 * to its default, builds `latest` (recomputeLatest), then hands off to
 * the left-pane sections and renderFilteredExplorer for the actual DOM
 * build.
 */
async function renderContigTable(records, binTablesByTool) {
  currentRecords = records;
  currentParams = defaultGlobalParams();
  currentMagFilters = window.ClannMAG.magFilters.defaultMagFilters();
  currentFilters = window.ClannMAG.filters.defaultFilters();
  allBinTablesByTool = binTablesByTool;
  activeTools = new Set(binTablesByTool ? binTablesByTool.keys() : []);

  await recomputeLatest(records, filterActiveBinTables());

  document.getElementById('empty').style.display = 'none';
  document.getElementById('explorer').style.display = 'flex';

  renderFiltersSection();
  renderMagFiltersSection();
  renderParamsSection();
  renderToolsSection();
  renderFilteredExplorer();
}

// Parsing runs in a Worker (src/workers/fasta-worker.js) rather than on the
// main thread — see docs/phase1-investigation.md "Performance follow-ups":
// same total work, but the page stays responsive and can show live
// progress instead of freezing for however long a large assembly takes.
function attachAuxiliaryData(records, coverageTable, krakenCalls) {
  if (coverageTable) {
    const depthsByContig = new Map(coverageTable.rows.map((r) => [r.contigId, r.depths]));
    for (const record of records) {
      const depths = depthsByContig.get(record.id);
      if (depths) record.coverageDepths = depths;
    }
  }
  if (krakenCalls) {
    const taxIdByContig = new Map(krakenCalls.filter((c) => c.classified).map((c) => [c.contigId, c.taxId]));
    for (const record of records) {
      const taxId = taxIdByContig.get(record.id);
      if (taxId != null) record.krakenTaxId = taxId;
    }
  }
}

// Marker-gene search is the slow part of loading an assembly (a hand-rolled
// seed-and-extend search per contig — see src/model/marker-genes.js), so it
// runs in its own pool of workers, one per available core, dispatched
// round-robin as fasta-worker.js streams contigs in. fasta-worker.js itself
// stays single-threaded (the streaming parse doesn't shard cleanly) but
// marker search is embarrassingly parallel across contigs.
const MARKER_POOL_SIZE = Math.max(1, Math.min(navigator.hardwareConcurrency || 4, 8));

/**
 * Union of every contigId referenced across all loaded bin tables (a
 * contig kept by any one tool's assignment counts), or null if no bin
 * table was loaded at all — in which case there's nothing to filter
 * against, so every contig in the FASTA is kept (the plain assembly-QC-
 * only workflow). Computed from the already-parsed tables (findAuxiliaryFiles
 * runs before loadAssembly is ever called — see initFilePicker), then
 * handed to fasta-worker.js so contigs absent from every bin table are
 * dropped during the stream itself rather than parsed, six-frame-
 * translated, and marker-searched only to be discarded afterward.
 */
function computeReferencedContigIds(binTablesByTool) {
  if (!binTablesByTool) return null;
  const ids = new Set();
  for (const assignments of binTablesByTool.values()) {
    for (const { contigId } of assignments) ids.add(contigId);
  }
  return ids;
}

function loadAssembly(file, binTablesByTool, coverageTable, krakenCalls) {
  const referencedContigIds = computeReferencedContigIds(binTablesByTool);
  return new Promise((resolve, reject) => {
    const worker = new Worker('src/workers/fasta-worker.js');
    const markerWorkers = Array.from({ length: MARKER_POOL_SIZE }, () => new Worker('src/workers/marker-search-worker.js'));
    const pendingById = new Map(); // requestId -> resolve
    const pendingPromises = [];
    let nextWorker = 0;
    let nextRequestId = 0;

    markerWorkers.forEach((w) => {
      w.onmessage = (e) => {
        const { id, markerHits } = e.data;
        const resolveOne = pendingById.get(id);
        if (resolveOne) { pendingById.delete(id); resolveOne(markerHits); }
      };
    });

    function dispatchMarkerSearch(record) {
      if (!record.frames) return;
      const id = nextRequestId++;
      const w = markerWorkers[nextWorker];
      nextWorker = (nextWorker + 1) % markerWorkers.length;
      const frames = record.frames;
      delete record.frames; // handed off; the worker gets its own copy via postMessage
      const done = new Promise((res) => pendingById.set(id, res)).then((markerHits) => {
        record.markerHits = markerHits;
      });
      pendingPromises.push(done);
      w.postMessage({ id, frames });
    }

    function terminateAll() {
      worker.terminate();
      markerWorkers.forEach((w) => w.terminate());
      setBusy(false);
    }

    const records = [];
    const t0 = performance.now();
    setBusy(true);
    showError(`Parsing ${file.name}… 0 contigs so far`);
    currentAssemblyFile = file;
    workingAssignmentInitialized = false;
    networkAlgorithm = 'ring';
    selectedMagId = null;
    selectedContigId = null;
    showUncontendedContigs = false;
    showConnectedMagContigs = false;
    decidedContigIds = new Set();

    worker.onmessage = async (e) => {
      const msg = e.data;
      if (msg.type === 'contig') {
        records.push(msg.record);
        dispatchMarkerSearch(msg.record);
      } else if (msg.type === 'progress') {
        showError(`Parsing ${file.name}… ${msg.contigsSoFar.toLocaleString()} contigs so far`);
      } else if (msg.type === 'done') {
        await Promise.all(pendingPromises);
        const elapsed = ((performance.now() - t0) / 1000).toFixed(1);
        showError(null);
        attachAuxiliaryData(records, coverageTable, krakenCalls);
        await renderContigTable(records, binTablesByTool);
        const excludedNote = msg.summary.excludedCount > 0
          ? ` · ${msg.summary.excludedCount.toLocaleString()} unreferenced contigs excluded`
          : '';
        document.getElementById('hMeta').textContent =
          `${msg.summary.contigCount.toLocaleString()} contigs · ${msg.summary.totalLength.toLocaleString()} bp · parsed in ${elapsed}s${excludedNote}`;
        document.getElementById('hTitle').textContent = file.name;
        terminateAll();
        resolve();
      } else if (msg.type === 'error') {
        showError(`Failed to parse ${file.name}: ${msg.message}`);
        terminateAll();
        reject(new Error(msg.message));
      }
    };
    worker.onerror = (err) => {
      showError(`Failed to parse ${file.name}: ${err.message}`);
      terminateAll();
      reject(err);
    };
    worker.postMessage({ file, referencedContigIds });
  });
}

/** Strips a filename down to its base name (no directories, no extension) for use as a tool label. */
function fileBaseName(name) {
  const withoutDir = name.split('/').pop();
  const dot = withoutDir.lastIndexOf('.');
  return dot > 0 ? withoutDir.slice(0, dot) : withoutDir;
}

/**
 * Among the files NOT identified as the assembly, content-sniff every one
 * and route it: any number of contig->bin tables (Phase 4/5), at most one
 * coverage table (brief's inputs list one coverage table, not per-tool),
 * and at most one Kraken2 per-contig call file (Phase 6). First match
 * wins for the single-instance inputs, matching how the assembly pick
 * itself already works; bin tables collect all matches, labelled by
 * filename since content alone can't name which tool produced a table.
 */
async function findAuxiliaryFiles(otherFiles) {
  const { sniff } = window.ClannMAG.sniff;
  const { parseContigBinTable } = window.ClannMAG.contigBinTable;
  const { parseCoverageTable } = window.ClannMAG.coverageTable;
  const { parseKraken2ContigCalls } = window.ClannMAG.kraken2Contigs;

  const binTablesByTool = new Map();
  let coverageTable = null;
  let krakenCalls = null;

  for (const file of otherFiles) {
    const text = await file.text();
    const format = sniff(text).format;
    if (format === 'contig-bin-table') {
      let label = fileBaseName(file.name);
      let suffix = 2;
      while (binTablesByTool.has(label)) label = `${fileBaseName(file.name)} (${suffix++})`;
      binTablesByTool.set(label, parseContigBinTable(text));
    } else if (format === 'coverage-table' && !coverageTable) {
      coverageTable = parseCoverageTable(text);
    } else if (format === 'kraken2-contigs' && !krakenCalls) {
      krakenCalls = parseKraken2ContigCalls(text);
    }
  }
  return { binTablesByTool: binTablesByTool.size ? binTablesByTool : null, coverageTable, krakenCalls };
}

// Fetched once and cached: build/03-taxonomy.js's output, used for the
// marker-gene taxonomic-consistency check (brief §Marker-gene
// identification module). ~420KB, small enough to fetch lazily on first
// use rather than upfront alongside the (much larger) marker search
// assets, which the Worker already loads in parallel with file parsing.
let taxonomyTreePromise = null;
function getTaxonomyTree() {
  if (!taxonomyTreePromise) {
    taxonomyTreePromise = window.ClannMAG.markerTaxonomy.loadTaxonomyTree('data/').catch((err) => {
      console.warn('marker-gene taxonomy tree failed to load:', err.message);
      return null;
    });
  }
  return taxonomyTreePromise;
}

/**
 * Combines every per-contig signal this app computes (brief §Outlier and
 * disagreement flagging: "overlaid with the cross-tool agreement signal
 * ... combined with the other signals") into one ranked list: composition/
 * coverage centroid distance (outliers.js), marker contribution
 * (unique/redundant, bin-summary.js), marker-gene taxonomic consistency
 * (marker-taxonomy.js, if the lineage tree loaded), Kraken2 disagreement
 * (bin-summary.js, if a per-contig call file was loaded), and cross-tool
 * agreement (bin-reconciliation.js, if 2+ bin tables were loaded) — each
 * "hot" signal increments a contig's flagCount, the primary sort key.
 *
 * The contig grouping used as "the bin" for centroid/marker-contribution
 * purposes is the reconciled putative MAG's full contig set (core +
 * disputed) when 2+ tools are loaded, or the single table's own bins
 * otherwise — the best available guess at "this genome" either way.
 */
function computeOutlierFlags(records, binTablesByTool, reconciliation, tree) {
  const recordsById = new Map(records.map((r) => [r.id, r]));
  const groups = [];

  if (reconciliation) {
    for (const mag of reconciliation.putativeMags) {
      groups.push({ groupLabel: mag.magId, contigIds: [...mag.coreContigIds, ...mag.disputedContigIds] });
    }
  } else if (binTablesByTool) {
    const [[, onlyTable]] = binTablesByTool;
    const byBin = new Map();
    for (const { contigId, binId } of onlyTable) {
      if (!byBin.has(binId)) byBin.set(binId, []);
      byBin.get(binId).push(contigId);
    }
    for (const [binId, contigIds] of byBin) groups.push({ groupLabel: binId, contigIds });
  }

  const { computeBinOutliers } = window.ClannMAG.outliers;
  const { computeMarkerContributions, computeKrakenDisagreement } = window.ClannMAG.binSummary;
  const { computeBinTaxonomicConsistency } = window.ClannMAG.markerTaxonomy;

  const agreementByContig = new Map();
  if (reconciliation) for (const c of reconciliation.contigAgreement) agreementByContig.set(c.contigId, c);

  const flags = [];
  for (const group of groups) {
    const groupContigs = group.contigIds.map((id) => recordsById.get(id)).filter(Boolean);
    if (groupContigs.length === 0) continue;

    const zScores = computeBinOutliers(groupContigs);
    const contributions = computeMarkerContributions(groupContigs);
    const krakenFlags = computeKrakenDisagreement(groupContigs);
    const taxConsistency = tree ? computeBinTaxonomicConsistency(groupContigs, tree) : null;

    for (const contig of groupContigs) {
      const z = zScores.get(contig.id);
      const contribution = contributions.get(contig.id);
      const taxDistance = (taxConsistency && taxConsistency.perContigDistance.get(contig.id)) ?? null;
      const krakenDisagrees = krakenFlags.get(contig.id) || false;
      const agreement = agreementByContig.get(contig.id) || null;

      // flagCount isn't computed here — it depends on the adjustable
      // z/taxonomic-distance thresholds (outliers.js's computeFlagCount),
      // so it's derived at render time from these raw values instead
      // (applyOutlierThresholds, below), letting a threshold change
      // re-rank the list without redoing this per-contig analysis.
      flags.push({
        contigId: contig.id, groupLabel: group.groupLabel,
        compositionZ: z.compositionZ, coverageZ: z.coverageZ, combinedZ: z.combinedZ,
        uniqueCount: contribution.uniqueFamilies.length, redundantCount: contribution.redundantFamilies.length,
        redundantOnly: contribution.redundantFamilies.length > 0 && contribution.uniqueFamilies.length === 0,
        taxDistance, krakenDisagrees, agreementFraction: agreement ? agreement.agreementFraction : null,
      });
    }
  }

  return flags;
}

/**
 * Applies the current (adjustable) outlier-flagging thresholds to raw
 * per-contig signals (computeOutlierFlags's output) and re-sorts — split
 * out so a threshold change in the Thresholds & parameters panel can
 * re-rank the outlier table without recomputing composition/coverage
 * z-scores, marker contributions, or taxonomic consistency from scratch.
 */
function applyOutlierThresholds(flags, outlierParams) {
  const { computeFlagCount } = window.ClannMAG.outliers;
  const withCounts = flags.map((f) => ({ ...f, flagCount: computeFlagCount(f, outlierParams) }));
  withCounts.sort((a, b) => b.flagCount - a.flagCount || b.combinedZ - a.combinedZ);
  return withCounts;
}

/**
 * Scoped, per the redesign, to whichever MAG is currently selected in the
 * picker table (the caller filters `flags` down to that MAG's network
 * neighborhood before calling this — see renderFilteredExplorer) rather
 * than showing every contig in the assembly at once.
 */
function renderOutlierCard(flags, { hasCoverage, hasTaxonomy, hasKraken, hasCrossTool }, selectedMagId) {
  if (!selectedMagId) {
    return `
      <div class="card">
        <h3>Outlier &amp; disagreement flagging</h3>
        <div class="hint">Select a putative MAG above to see outlier/disagreement flags for its contigs.</div>
      </div>
    `;
  }
  const ROW_LIMIT = 150;
  const rows = flags
    .slice(0, ROW_LIMIT)
    .map((f) => `<tr>
      <td>${f.contigId}</td>
      <td>${f.groupLabel}</td>
      <td class="num">${f.compositionZ.toFixed(2)}</td>
      <td class="num">${f.coverageZ === null ? '<span class="hint">n/a</span>' : f.coverageZ.toFixed(2)}</td>
      <td class="num">${f.uniqueCount}</td>
      <td class="num">${f.redundantCount}</td>
      <td class="num">${f.taxDistance === null ? '<span class="hint">n/a</span>' : f.taxDistance}</td>
      <td>${f.krakenDisagrees ? '⚠' : ''}</td>
      <td class="num">${f.agreementFraction === null ? '<span class="hint">n/a</span>' : `${(f.agreementFraction * 100).toFixed(0)}%`}</td>
      <td class="num"><strong>${f.flagCount}</strong></td>
    </tr>`)
    .join('');

  const notes = [
    'composition Z is always available',
    hasCoverage ? 'coverage Z from the loaded coverage table' : 'coverage Z: n/a, no coverage table loaded',
    hasTaxonomy ? 'marker taxonomic distance from the loaded lineage table' : 'marker taxonomic distance: n/a, lineage table failed to load',
    hasKraken ? 'Kraken2 disagreement from the loaded per-contig calls' : 'Kraken2 disagreement: n/a, no per-contig Kraken2 file loaded',
    hasCrossTool ? 'cross-tool agreement from the reconciled bin tables' : 'cross-tool agreement: n/a, load 2+ bin tables to enable',
  ];

  return `
    <div class="card">
      <h3>Outlier &amp; disagreement flagging — ${selectedMagId}'s neighborhood</h3>
      <div class="row-count">${flags.length.toLocaleString()} contigs scored${flags.length > ROW_LIMIT ? `, showing the top ${ROW_LIMIT} by flag count` : ''} &middot; ${notes.join(' &middot; ')}</div>
      <div class="table-wrap scroll-panel">
        <table class="data-table">
          <thead><tr>
            <th>Contig</th><th>Bin / MAG</th>
            <th title="Standard deviations from this bin's composition centroid">Comp Z</th>
            <th title="Standard deviations from this bin's coverage centroid">Cov Z</th>
            <th title="Marker families found only on this contig within its bin">Unique</th>
            <th title="Marker families also found elsewhere in this bin">Redundant</th>
            <th title="How far this contig's marker provenance sits from the rest of its bin's consensus lineage">Tax dist.</th>
            <th title="This contig's Kraken2 call disagrees with its bin's majority call">Kraken</th>
            <th title="Cross-tool agreement fraction (Phase 5)">Agreement</th>
            <th title="Count of signals flagged for this contig">Flags</th>
          </tr></thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
    </div>
  `;
}

function renderContigMatchNotices(notices) {
  const el = document.getElementById('contigMatchNotices');
  if (!el) return;
  el.innerHTML = notices
    .map((n) => `<div class="hint${n.level === 'warning' ? ' hint-warn' : ''}">${n.level === 'warning' ? '⚠' : '✓'} ${n.message}</div>`)
    .join('');
}

/**
 * Best-attempt contig ID matching (src/model/contig-id-matching.js) run
 * once per load, before the real streaming parse: some binning tools
 * (CONCOCT's default contig-splitting step, most concretely) rewrite
 * contig IDs before clustering, which otherwise makes every one of that
 * tool's rows silently fail to match the assembly and contribute nothing
 * to bin summaries or reconciliation — not because the tool was run
 * against a different assembly, but because of a known, mechanical naming
 * convention. A quick header-only pre-scan of the assembly (scanContigIds
 * — no per-contig stats, unlike the real parse) gets the real contig ID
 * set needed to detect and fix this per tool; a tool whose IDs still
 * don't resolve well after checking known conventions gets a visible
 * warning instead, since that's the genuine "this might be a different
 * assembly version" signal the brief's provenance concern is about, and
 * guessing further would just produce confident-looking wrong results.
 */
async function bestAttemptMatchBinTables(assemblyFile, binTablesByTool) {
  const { scanContigIds } = window.ClannMAG.fastaIndex;
  const { bestAttemptRemapAssignments, HIGH_MATCH_THRESHOLD } = window.ClannMAG.contigIdMatching;

  setBusy(true);
  showError('Checking contig IDs against the assembly…');
  let referenceIds;
  try {
    referenceIds = await scanContigIds(assemblyFile);
  } finally {
    setBusy(false);
    showError(null);
  }

  const remapped = new Map();
  const notices = [];
  for (const [tool, assignments] of binTablesByTool) {
    const { assignments: fixedAssignments, report } = bestAttemptRemapAssignments(assignments, referenceIds);
    remapped.set(tool, fixedAssignments);
    if (report.applied) {
      notices.push({
        tool, level: 'fixed',
        message: `${tool}: table used ${report.patternLabel}; corrected automatically `
          + `(${(report.matchRateBefore * 100).toFixed(0)}% → ${(report.matchRateAfter * 100).toFixed(0)}% matched the assembly)`
          + `${report.collapsedCount ? `, ${report.collapsedCount.toLocaleString()} split contig(s) merged back via majority vote across their parts` : ''}.`,
      });
    } else if (report.matchRateBefore < HIGH_MATCH_THRESHOLD) {
      notices.push({
        tool, level: 'warning',
        message: `${tool}: only ${(report.matchRateAfter * 100).toFixed(0)}% of its contig IDs match the loaded assembly, `
          + `even after checking known naming conventions — it may have been run against a different assembly version. `
          + `Its bin calls will mostly be skipped as unmatched; consider excluding it from the comparison.`,
      });
    }
  }
  renderContigMatchNotices(notices);
  return remapped;
}

function initFilePicker() {
  const input = document.getElementById('folder-input');
  const openButtons = [document.getElementById('uploadBtn'), document.getElementById('emptyOpen')];
  openButtons.forEach((btn) => btn && btn.addEventListener('click', () => input.click()));
  input.addEventListener('change', async () => {
    const files = [...input.files];
    if (files.length === 0) return;
    renderContigMatchNotices([]);

    let assemblyFile = null;
    for (const file of files) {
      if (await looksLikeFasta(file)) { assemblyFile = file; break; }
    }
    if (!assemblyFile) {
      showError('No file among your selection looks like a FASTA assembly (checked for a leading ">" or gzip magic bytes).');
      return;
    }
    const otherFiles = files.filter((f) => f !== assemblyFile);
    const { binTablesByTool, coverageTable, krakenCalls } = otherFiles.length
      ? await findAuxiliaryFiles(otherFiles)
      : { binTablesByTool: null, coverageTable: null, krakenCalls: null };

    const matchedBinTablesByTool = binTablesByTool
      ? await bestAttemptMatchBinTables(assemblyFile, binTablesByTool)
      : null;

    try {
      await loadAssembly(assemblyFile, matchedBinTablesByTool, coverageTable, krakenCalls);
    } catch {
      // already surfaced via showError inside loadAssembly
    }
  });
}

/**
 * Click-to-sort for every `table.data-table` in the page (contig table,
 * bin summaries, comparison/outlier views, etc. — all of them already
 * share this one class). Wired via a single delegated listener on
 * `document` rather than per-table, since tables here are re-rendered
 * wholesale (innerHTML swaps) as data changes — a delegated listener
 * keeps working across re-renders with nothing to re-attach. The CSS for
 * the sorted-column arrow (`.sorted-asc`/`.sorted-desc`) already existed
 * (styles/main.css) with no JS behind it; this is that missing half.
 *
 * First click on a header sorts ascending, a second click on the same
 * header flips to descending, and clicking a different header restarts
 * at ascending for that column (clearing the previous column's arrow).
 * Column type (numeric vs text) is inferred per sort from the cells
 * actually present, so it stays correct as rows are added/removed/edited
 * between clicks rather than being decided once up front.
 */
function initSortableTables() {
  document.addEventListener('click', (e) => {
    const th = e.target.closest('table.data-table thead th');
    if (!th) return;
    const table = th.closest('table');
    const tbody = table.querySelector('tbody');
    if (!tbody) return;
    const headerRow = th.parentElement;
    const headerCells = [...headerRow.children];
    const colIndex = headerCells.indexOf(th);
    if (colIndex < 0) return;

    const nextDir = th.classList.contains('sorted-asc') ? 'desc' : 'asc';
    headerCells.forEach((cell) => cell.classList.remove('sorted-asc', 'sorted-desc'));
    th.classList.add(nextDir === 'asc' ? 'sorted-asc' : 'sorted-desc');

    const rows = [...tbody.rows];
    const cellValue = (row) => {
      const text = (row.cells[colIndex] ? row.cells[colIndex].textContent : '').trim();
      const num = text === '' ? NaN : Number(text.replace(/[,%]/g, ''));
      return { text, num };
    };
    const values = rows.map(cellValue);
    const numericCount = values.filter((v) => !Number.isNaN(v.num)).length;
    const isNumeric = numericCount >= values.length * 0.5; // majority-numeric column, e.g. a few "n/a" cells mixed in with numbers

    const withValues = rows.map((row, i) => ({ row, value: values[i] }));
    withValues.sort((a, b) => {
      const aEmpty = isNumeric ? Number.isNaN(a.value.num) : a.value.text === '';
      const bEmpty = isNumeric ? Number.isNaN(b.value.num) : b.value.text === '';
      if (aEmpty || bEmpty) return aEmpty - bEmpty; // blank/n-a cells always sink to the bottom, regardless of direction
      const cmp = isNumeric ? a.value.num - b.value.num : a.value.text.localeCompare(b.value.text);
      return nextDir === 'asc' ? cmp : -cmp;
    });
    for (const { row } of withValues) tbody.appendChild(row);
  });
}

/**
 * Delegated (survives renderReconciliationCard's wholesale innerHTML
 * rebuilds — same reasoning as initSortableTables above) click handler for
 * the MAG picker table's row-select buttons. Clicking the already-selected
 * MAG deselects it (toggle), clearing the network/evidence panel back to
 * the "pick a MAG" hint state.
 */
function initMagPicker() {
  document.addEventListener('click', (e) => {
    const btn = e.target.closest('.mag-picker-select');
    if (!btn || !latest) return;
    const magId = btn.dataset.magId;
    selectedMagId = selectedMagId === magId ? null : magId;
    selectedContigId = null;
    renderFilteredExplorer();
  });
}

document.addEventListener('DOMContentLoaded', () => {
  initTheme();
  initFilePicker();
  initSortableTables();
  initMagPicker();
});
})();
