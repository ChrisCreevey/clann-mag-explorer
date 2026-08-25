(function () {
  'use strict';

// App shell. Loads an assembly FASTA (streamed through fasta-index.js in a
// Worker, with marker-gene search) and, per Phase 4/5, any number of
// contig->bin tables alongside it — content-sniffed from whatever else was
// selected, not by filename (the tool label per table IS taken from its
// filename, since content alone can't name which binning tool produced
// it). One table gets a single per-bin summary (Phase 4); two or more
// also get cross-tool reconciliation (Phase 5): bins matched across tools
// by contig overlap, core/disputed contig sets, and a side-by-side view.
// Phase 6 adds a combined outlier/disagreement view; Phase 7 (this
// module's newest section, see renderInteractiveSection below) adds a
// live-editable working bin assignment — a scatter plot with rectangular
// drag-select (brief's simplified lasso, see scatter-geometry.js), move/
// merge/split-via-reassignment actions, and bin summaries that
// recalculate immediately off the current in-session edits rather than
// the original loaded tables. Filters/search land in a later phase.

const THEME_KEY = 'clann-mag-explorer-theme';

// Session-level interactive-reassignment state — deliberately module-level
// (not passed around as parameters) since it needs to survive across
// scatter-plot re-renders and outlives any single render() call. Reset on
// every fresh assembly load (loadAssembly), including the unsaved-work
// flag, matching the brief's "the tool must intercept navigation away
// ... once any reassignment has been made in the session" — a new session
// starts clean.
let workingAssignment = new Map();
let hasUnsavedReassignments = false;
let scatterPlotHandle = null;

// Phase 9 export: the loaded assembly's own File (re-sliced via
// Blob.slice for per-MAG FASTA extraction, never re-uploaded) and its
// parsed records, kept at module scope alongside workingAssignment since
// export reads from the same live session state. Reset on every fresh
// load, same as workingAssignment above — a reload loses the File
// reference regardless (brief's accepted, documented limitation; see
// docs/phase1-investigation.md Phase 1), so there's nothing to persist.
let currentAssemblyFile = null;
let currentRecords = [];

window.addEventListener('beforeunload', (e) => {
  if (!hasUnsavedReassignments) return;
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

function renderBinSummaryCard(records, binAssignments, toolLabel) {
  const { computeBinSummaries } = window.ClannMAG.binSummary;
  const { summaries, unmatchedContigIds } = computeBinSummaries(records, binAssignments);

  const rows = summaries
    .map((b) => `<tr>
      <td>${b.binId}</td>
      <td class="num">${b.contigCount.toLocaleString()}</td>
      <td class="num">${b.totalLength.toLocaleString()}</td>
      <td class="num">${b.n50.toLocaleString()}</td>
      <td class="num">${b.l50.toLocaleString()}</td>
      <td class="num">${(b.meanGc * 100).toFixed(1)}%</td>
      <td class="num">${b.completeness.toFixed(1)}%</td>
      <td class="num">${b.redundancy.toFixed(1)}%</td>
      <td>${formatMimagTier(b.mimagTier)}</td>
    </tr>`)
    .join('');

  const unmatchedNote = unmatchedContigIds.length
    ? `<div class="hint">${unmatchedContigIds.length.toLocaleString()} contig ID(s) in the bin table were not found in the loaded assembly and were skipped (e.g. ${unmatchedContigIds.slice(0, 3).join(', ')}${unmatchedContigIds.length > 3 ? ', …' : ''}).</div>`
    : '';

  return `
    <div class="card">
      <h3>Bin summaries${toolLabel ? ` — ${toolLabel}` : ''}</h3>
      <div class="row-count">${summaries.length.toLocaleString()} bins, largest first &middot; completeness/redundancy from the built-in marker-gene search (40 families) &middot; MIMAG-style tier is a completeness/contamination proxy only, not the full standard (no rRNA/tRNA check)</div>
      <div class="table-wrap scroll-panel">
        <table class="data-table">
          <thead><tr>
            <th>Bin</th><th>Contigs</th><th>Length</th><th>N50</th><th>L50</th><th>Mean GC</th>
            <th title="Fraction of the 40 marker families found anywhere in this bin">Completeness</th>
            <th title="Fraction of found families that appear on more than one contig — a proxy for contamination">Redundancy</th>
            <th>Tier</th>
          </tr></thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
      ${unmatchedNote}
    </div>
  `;
}

/**
 * Phase 5: matches bins across two or more loaded tools by contig
 * overlap and renders the reconciled view — putative MAGs (side by side
 * across tools), and the ranked disputed-contig list.
 */
function renderReconciliationCard(records, result) {
  const { computeCompletenessRedundancy, mimagTier } = window.ClannMAG.binSummary;
  const recordsById = new Map(records.map((r) => [r.id, r]));
  const tools = result.tools;

  const magRows = result.putativeMags
    .map((mag) => {
      const coreRecords = mag.coreContigIds.map((id) => recordsById.get(id)).filter(Boolean);
      const { completeness, redundancy } = computeCompletenessRedundancy(coreRecords);
      const tier = mimagTier(completeness, redundancy);
      const toolCells = tools
        .map((tool) => {
          const member = mag.members.find((m) => m.tool === tool);
          return member
            ? `<td>${member.binId} <span class="hint">(${member.contigCount.toLocaleString()})</span></td>`
            : '<td class="hint">—</td>';
        })
        .join('');
      return `<tr>
        <td>${mag.magId}</td>
        <td class="num">${mag.coreContigIds.length.toLocaleString()}</td>
        <td class="num">${mag.disputedContigIds.length.toLocaleString()}</td>
        <td class="num">${completeness.toFixed(1)}%</td>
        <td class="num">${redundancy.toFixed(1)}%</td>
        <td>${formatMimagTier(tier)}</td>
        ${toolCells}
      </tr>`;
    })
    .join('');

  const DISPUTED_ROW_LIMIT = 200;
  const disputedRows = result.disputedContigsRanked
    .slice(0, DISPUTED_ROW_LIMIT)
    .map((c) => {
      const votesStr = tools.map((t) => `${t}: ${c.votes[t] ?? '<span class="hint">unbinned</span>'}`).join(' &middot; ');
      return `<tr>
        <td>${c.contigId}</td>
        <td class="num">${(c.agreementFraction * 100).toFixed(0)}%</td>
        <td class="num">${c.totalVotes}</td>
        <td>${votesStr}</td>
      </tr>`;
    })
    .join('');

  return `
    <div class="card">
      <h3>Cross-tool reconciliation</h3>
      <div class="row-count">${tools.length} tools loaded (${tools.join(', ')}) &middot; ${result.putativeMags.length.toLocaleString()} putative MAGs matched by contig overlap (reciprocal best hit, min Jaccard 0.1) &middot; completeness/redundancy below are computed from each MAG's high-confidence core (unanimous-agreement) contigs only</div>
      <div class="table-wrap scroll-panel">
        <table class="data-table">
          <thead><tr>
            <th>Putative MAG</th>
            <th title="Contigs every voting tool agrees belong to this MAG">Core</th>
            <th title="Contigs assigned here by some but not all voting tools">Disputed</th>
            <th>Completeness</th>
            <th>Redundancy</th>
            <th>Tier</th>
            ${tools.map((t) => `<th>${t}</th>`).join('')}
          </tr></thead>
          <tbody>${magRows}</tbody>
        </table>
      </div>
      <h4>Disputed contigs, most split first</h4>
      <div class="row-count">${result.disputedContigsRanked.length.toLocaleString()} contigs where loaded tools disagree${result.disputedContigsRanked.length > DISPUTED_ROW_LIMIT ? `, showing the first ${DISPUTED_ROW_LIMIT}` : ''}</div>
      <div class="table-wrap scroll-panel">
        <table class="data-table">
          <thead><tr><th>Contig</th><th>Agreement</th><th>Tools voting</th><th>Votes by tool</th></tr></thead>
          <tbody>${disputedRows}</tbody>
        </table>
      </div>
    </div>
  `;
}

/**
 * Phase 7: builds and wires the live-editable reassignment section — a
 * scatter plot (brief's "composition, coverage, GC, length" axes; here
 * per-contig scalars, not the raw 136-dim composition vector, since a
 * literal scatter can't plot that many dimensions at once — Phase 6's
 * compositionZ is available as a derived axis instead), rectangular
 * drag-select, move/merge/new-bin actions, and a working-bins summary
 * table that recalculates immediately from bin-summary.js against the
 * current in-session assignment (working-assignment.js), not the
 * original loaded tables.
 */
function initInteractiveSection(records, binTablesByTool, reconciliationResult) {
  const {
    deriveInitialAssignment, assignmentToRows, reassignContigs, generateNewBinId, listBinIds,
  } = window.ClannMAG.workingAssignment;
  const { computeBinSummaries } = window.ClannMAG.binSummary;
  const { createScatterPlot } = window.ClannMAG.scatter;

  workingAssignment = deriveInitialAssignment(binTablesByTool, reconciliationResult);
  hasUnsavedReassignments = false;

  const AXES = {
    gcContent: { label: 'GC%', get: (r) => r.gcContent * 100 },
    length: { label: 'Length (log₁₀ bp)', get: (r) => Math.log10(Math.max(1, r.length)) },
    gcSkew: { label: 'GC skew', get: (r) => r.gcSkew },
    codingDensity: { label: 'Coding density %', get: (r) => r.codingDensity * 100 },
    coverage: {
      label: 'Mean coverage depth',
      get: (r) => (r.coverageDepths ? r.coverageDepths.reduce((a, b) => a + b, 0) / r.coverageDepths.length : 0),
    },
  };
  const hasCoverage = records.some((r) => r.coverageDepths);
  const axisKeys = Object.keys(AXES).filter((k) => k !== 'coverage' || hasCoverage);
  const axisOptionsHtml = (selectedKey) =>
    axisKeys.map((k) => `<option value="${k}" ${k === selectedKey ? 'selected' : ''}>${AXES[k].label}</option>`).join('');

  const card = document.getElementById('interactive-card');
  card.innerHTML = `
    <h3>Refine bins (interactive)</h3>
    <div class="row-count">Drag on the plot to select a cluster of contigs (a simplified rectangular lasso), then move them to a bin below. Reassignments recalculate bin stats immediately; export of the revised assignment table is a later phase.</div>
    <div class="row"><label>X axis</label><select id="scatterX">${axisOptionsHtml('gcContent')}</select></div>
    <div class="row"><label>Y axis</label><select id="scatterY">${axisOptionsHtml('length')}</select></div>
    <div id="scatterContainer"></div>
    <div class="row" id="selectionRow" style="display:none"><label id="selectionLabel"></label></div>
    <div class="row" id="actionRow" style="display:none">
      <select id="targetBinSelect"></select>
      <button class="act" id="applyMoveBtn" type="button">Move selected</button>
      <button class="act" id="clearSelectionBtn" type="button">Clear selection</button>
    </div>
    <h4>Current working bins</h4>
    <div class="row-count" id="workingBinsNote"></div>
    <div class="table-wrap scroll-panel" id="workingBinsWrap"></div>
  `;

  function buildDataPoints(xKey, yKey) {
    return records.map((r) => ({
      id: r.id, x: AXES[xKey].get(r), y: AXES[yKey].get(r),
      colorKey: workingAssignment.get(r.id) || 'unbinned',
    }));
  }

  function renderWorkingBins() {
    const { summaries } = computeBinSummaries(records, assignmentToRows(workingAssignment));
    const rows = summaries
      .map((b) => `<tr>
        <td>${b.binId}</td>
        <td class="num">${b.contigCount.toLocaleString()}</td>
        <td class="num">${b.totalLength.toLocaleString()}</td>
        <td class="num">${b.completeness.toFixed(1)}%</td>
        <td class="num">${b.redundancy.toFixed(1)}%</td>
        <td>${formatMimagTier(b.mimagTier)}</td>
      </tr>`)
      .join('');
    document.getElementById('workingBinsWrap').innerHTML = `
      <table class="data-table">
        <thead><tr><th>Bin</th><th>Contigs</th><th>Length</th><th>Completeness</th><th>Redundancy</th><th>Tier</th></tr></thead>
        <tbody>${rows}</tbody>
      </table>
    `;
    const unassignedCount = records.length - workingAssignment.size;
    document.getElementById('workingBinsNote').textContent =
      `${summaries.length.toLocaleString()} working bins` +
      (unassignedCount > 0 ? ` · ${unassignedCount.toLocaleString()} contig(s) not yet assigned to any bin` : '');
  }

  function refreshTargetBinOptions() {
    const select = document.getElementById('targetBinSelect');
    const previousValue = select.value;
    select.innerHTML = listBinIds(workingAssignment).map((b) => `<option value="${b}">${b}</option>`).join('') +
      '<option value="__new__">+ New bin…</option>';
    if ([...select.options].some((o) => o.value === previousValue)) select.value = previousValue;
  }

  function updateSelectionUI(selectedIds) {
    const selectionRow = document.getElementById('selectionRow');
    const actionRow = document.getElementById('actionRow');
    if (selectedIds.length === 0) {
      selectionRow.style.display = 'none';
      actionRow.style.display = 'none';
      return;
    }
    selectionRow.style.display = '';
    actionRow.style.display = '';
    document.getElementById('selectionLabel').textContent = `${selectedIds.length.toLocaleString()} contig(s) selected`;
    refreshTargetBinOptions();
  }

  scatterPlotHandle = createScatterPlot(
    document.getElementById('scatterContainer'),
    buildDataPoints('gcContent', 'length'),
    { onSelectionChange: updateSelectionUI }
  );

  function rerenderScatterColors() {
    const xKey = document.getElementById('scatterX').value;
    const yKey = document.getElementById('scatterY').value;
    scatterPlotHandle.setDataPoints(buildDataPoints(xKey, yKey));
  }
  document.getElementById('scatterX').addEventListener('change', rerenderScatterColors);
  document.getElementById('scatterY').addEventListener('change', rerenderScatterColors);

  document.getElementById('applyMoveBtn').addEventListener('click', () => {
    const selectedIds = scatterPlotHandle.getSelection();
    if (selectedIds.length === 0) return;
    let targetBinId = document.getElementById('targetBinSelect').value;
    if (targetBinId === '__new__') {
      const name = window.prompt('New bin name:', generateNewBinId(workingAssignment));
      if (!name) return;
      targetBinId = name;
    }
    workingAssignment = reassignContigs(workingAssignment, selectedIds, targetBinId);
    hasUnsavedReassignments = true;
    rerenderScatterColors();
    scatterPlotHandle.clearSelection();
    renderWorkingBins();
  });

  document.getElementById('clearSelectionBtn').addEventListener('click', () => {
    scatterPlotHandle.clearSelection();
  });

  renderWorkingBins();
}

/**
 * Phase 8: comparison views and QC across the full reconciled set (brief
 * §Comparison and QC) — good-vs-bad side-by-side contig scatter,
 * completeness/contamination scatter across every putative MAG coloured
 * by supporting tools, and the pairwise redundancy check flagging likely
 * duplicate genomes. Only meaningful once cross-tool reconciliation has
 * run (needs `putativeMags`), so callers only invoke this when more than
 * one bin table was loaded.
 */
function initQcSection(records, reconciliationResult) {
  const { computeBinSummaries } = window.ClannMAG.binSummary;
  const { pickComparisonBins, buildMagQcPoints } = window.ClannMAG.qcComparison;
  const { computeMagRedundancy } = window.ClannMAG.magRedundancy;
  const { createScatterPlot } = window.ClannMAG.scatter;

  const { putativeMags } = reconciliationResult;
  const magAssignmentRows = putativeMags.flatMap((mag) =>
    [...mag.coreContigIds, ...mag.disputedContigIds].map((contigId) => ({ contigId, binId: mag.magId }))
  );
  const { summaries: magSummaries } = computeBinSummaries(records, magAssignmentRows);

  const card = document.getElementById('qc-card');
  card.innerHTML = `
    <h3>Comparison and QC across putative MAGs</h3>
    <h4>Good bin vs. bad bin</h4>
    <div class="row-count" id="qcComparisonNote"></div>
    <div class="qc-comparison-row">
      <div><div class="row-count" id="qcGoodLabel"></div><div id="qcGoodScatter"></div></div>
      <div><div class="row-count" id="qcBadLabel"></div><div id="qcBadScatter"></div></div>
    </div>
    <h4>Completeness vs. contamination, all putative MAGs</h4>
    <div class="row-count">${magSummaries.length.toLocaleString()} putative MAG(s), coloured by which tool(s) support each one. X axis: completeness %. Y axis: redundancy % (contamination proxy).</div>
    <div id="magQcScatter"></div>
    <h4>Possible duplicate genomes</h4>
    <div class="row-count" id="redundancyNote"></div>
    <div class="table-wrap scroll-panel" id="redundancyWrap"></div>
  `;

  const { good, bad } = pickComparisonBins(magSummaries);
  document.getElementById('qcComparisonNote').textContent = good && bad
    ? 'The MAG with the best completeness-minus-redundancy score, plotted against the worst, for a direct visual contrast.'
    : 'Need at least two putative MAGs with 2+ contigs each to show a comparison.';

  function contigScatterPoints(magId) {
    const contigIds = new Set(magAssignmentRows.filter((r) => r.binId === magId).map((r) => r.contigId));
    return records
      .filter((r) => contigIds.has(r.id))
      .map((r) => ({ id: r.id, x: r.gcContent * 100, y: Math.log10(Math.max(1, r.length)), colorKey: magId }));
  }

  if (good && bad) {
    document.getElementById('qcGoodLabel').textContent =
      `${good.binId} — completeness ${good.completeness.toFixed(1)}%, redundancy ${good.redundancy.toFixed(1)}%`;
    document.getElementById('qcBadLabel').textContent =
      `${bad.binId} — completeness ${bad.completeness.toFixed(1)}%, redundancy ${bad.redundancy.toFixed(1)}%`;
    createScatterPlot(document.getElementById('qcGoodScatter'), contigScatterPoints(good.binId),
      { width: 320, height: 260, onSelectionChange: () => {} });
    createScatterPlot(document.getElementById('qcBadScatter'), contigScatterPoints(bad.binId),
      { width: 320, height: 260, onSelectionChange: () => {} });
  }

  createScatterPlot(document.getElementById('magQcScatter'), buildMagQcPoints(magSummaries, putativeMags),
    { width: 640, height: 360, onSelectionChange: () => {} });

  const redundancyPairs = computeMagRedundancy(records, putativeMags);
  const flaggedPairs = redundancyPairs.filter((p) => p.likelyDuplicate);
  document.getElementById('redundancyNote').textContent = flaggedPairs.length > 0
    ? `${flaggedPairs.length.toLocaleString()} pair(s) of putative MAGs have near-identical composition — worth checking whether they're the same organism split across bins.`
    : 'No putative MAG pairs look like the same organism split across separate bins (by composition similarity).';
  document.getElementById('redundancyWrap').innerHTML = redundancyPairs.length === 0 ? '' : `
    <table class="data-table">
      <thead><tr><th>MAG A</th><th>MAG B</th><th>Composition similarity</th><th>Mean GC difference</th><th>Flag</th></tr></thead>
      <tbody>${redundancyPairs.map((p) => `<tr>
        <td>${p.magIdA}</td>
        <td>${p.magIdB}</td>
        <td class="num">${p.compositionSimilarity.toFixed(3)}</td>
        <td class="num">${(p.meanGcDiff * 100).toFixed(2)}%</td>
        <td>${p.likelyDuplicate ? '⚠ possible duplicate' : ''}</td>
      </tr>`).join('')}</tbody>
    </table>
  `;
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
    const { summaries } = computeBinSummaries(currentRecords, assignmentToRows(workingAssignment));
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

async function renderContigTable(records, binTablesByTool) {
  currentRecords = records;
  const sorted = [...records].sort((a, b) => b.length - a.length);
  const totalLength = sorted.reduce((sum, r) => sum + r.length, 0);
  const n50 = computeN50(sorted.map((r) => r.length), totalLength);
  const meanGc = sorted.length ? sorted.reduce((sum, r) => sum + r.gcContent, 0) / sorted.length : 0;
  const contigsWithMarkers = sorted.filter((r) => r.markerHits && r.markerHits.length > 0).length;
  const distinctFamiliesHit = new Set(sorted.flatMap((r) => (r.markerHits || []).map((h) => h.family))).size;

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

  const tools = binTablesByTool ? [...binTablesByTool.keys()] : [];
  const perToolBinCards = tools
    .map((tool) => renderBinSummaryCard(records, binTablesByTool.get(tool), tools.length > 1 ? tool : null))
    .join('');

  let reconciliationResult = null;
  if (tools.length > 1) {
    reconciliationResult = window.ClannMAG.binReconciliation.reconcileBins(binTablesByTool);
  }
  const reconciliationCard = reconciliationResult ? renderReconciliationCard(records, reconciliationResult) : '';

  let outlierCard = '';
  if (tools.length > 0) {
    const tree = await getTaxonomyTree();
    const flags = computeOutlierFlags(records, binTablesByTool, reconciliationResult, tree);
    outlierCard = renderOutlierCard(flags, {
      hasCoverage: records.some((r) => r.coverageDepths),
      hasTaxonomy: tree !== null,
      hasKraken: records.some((r) => r.krakenTaxId != null),
      hasCrossTool: reconciliationResult !== null,
    });
  }

  explorer.innerHTML = `
    <div class="card">
      <h3>Assembly summary</h3>
      <div class="row"><label>Contigs</label><strong>${sorted.length.toLocaleString()}</strong></div>
      <div class="row"><label>Total length</label><strong>${totalLength.toLocaleString()} bp</strong></div>
      <div class="row"><label>N50</label><strong>${n50.toLocaleString()} bp</strong></div>
      <div class="row"><label>Mean GC</label><strong>${(meanGc * 100).toFixed(1)}%</strong></div>
      <div class="row"><label>Contigs with marker genes</label><strong>${contigsWithMarkers.toLocaleString()}</strong></div>
      <div class="row"><label>Distinct marker families hit</label><strong>${distinctFamiliesHit} / 40</strong></div>
    </div>
    ${reconciliationCard}
    ${outlierCard}
    ${tools.length > 0 ? '<div class="card" id="interactive-card"></div>' : ''}
    ${reconciliationResult ? '<div class="card" id="qc-card"></div>' : ''}
    ${tools.length > 0 ? '<div class="card" id="export-card"></div>' : ''}
    ${perToolBinCards}
    <div class="card">
      <h3>Per-contig properties</h3>
      <div class="row-count">${sorted.length.toLocaleString()} contigs, longest first</div>
      <div class="table-wrap scroll-panel">
        <table class="data-table">
          <thead><tr><th>Contig</th><th>Length</th><th>GC%</th><th>GC skew</th><th>Coding density</th><th>Marker genes</th><th title="Non-uniform FASTA line wrapping">⚠</th></tr></thead>
          <tbody>${rows}</tbody>
        </table>
      </div>
    </div>
  `;
  document.getElementById('empty').style.display = 'none';
  explorer.style.display = 'flex';

  if (tools.length > 0) initInteractiveSection(records, binTablesByTool, reconciliationResult);
  if (reconciliationResult) initQcSection(records, reconciliationResult);
  if (tools.length > 0) initExportSection();
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

      let flagCount = 0;
      if (z.combinedZ > 2) flagCount++;
      if (contribution.redundantFamilies.length > 0 && contribution.uniqueFamilies.length === 0) flagCount++;
      if (taxDistance !== null && taxDistance > 0) flagCount++;
      if (krakenDisagrees) flagCount++;
      if (agreement && agreement.agreementFraction < 1) flagCount++;

      flags.push({
        contigId: contig.id, groupLabel: group.groupLabel,
        compositionZ: z.compositionZ, coverageZ: z.coverageZ, combinedZ: z.combinedZ,
        uniqueCount: contribution.uniqueFamilies.length, redundantCount: contribution.redundantFamilies.length,
        taxDistance, krakenDisagrees, agreementFraction: agreement ? agreement.agreementFraction : null,
        flagCount,
      });
    }
  }

  flags.sort((a, b) => b.flagCount - a.flagCount || b.combinedZ - a.combinedZ);
  return flags;
}

function renderOutlierCard(flags, { hasCoverage, hasTaxonomy, hasKraken, hasCrossTool }) {
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
      <h3>Outlier &amp; disagreement flagging</h3>
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

function initFilePicker() {
  const input = document.getElementById('folder-input');
  const openButtons = [document.getElementById('uploadBtn'), document.getElementById('emptyOpen')];
  openButtons.forEach((btn) => btn && btn.addEventListener('click', () => input.click()));
  input.addEventListener('change', async () => {
    const files = [...input.files];
    if (files.length === 0) return;

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

    try {
      await loadAssembly(assemblyFile, binTablesByTool, coverageTable, krakenCalls);
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

document.addEventListener('DOMContentLoaded', () => {
  initTheme();
  initFilePicker();
  initSortableTables();
});
})();
