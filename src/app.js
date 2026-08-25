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
// Filters and interactive reassignment land in later phases.

const THEME_KEY = 'clann-mag-explorer-theme';

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
function renderReconciliationCard(records, binTablesByTool) {
  const { reconcileBins } = window.ClannMAG.binReconciliation;
  const { computeCompletenessRedundancy, mimagTier } = window.ClannMAG.binSummary;
  const recordsById = new Map(records.map((r) => [r.id, r]));
  const result = reconcileBins(binTablesByTool);
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

function renderContigTable(records, binTablesByTool) {
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
  const reconciliationCard = tools.length > 1 ? renderReconciliationCard(records, binTablesByTool) : '';

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
}

// Parsing runs in a Worker (src/workers/fasta-worker.js) rather than on the
// main thread — see docs/phase1-investigation.md "Performance follow-ups":
// same total work, but the page stays responsive and can show live
// progress instead of freezing for however long a large assembly takes.
function loadAssembly(file, binTablesByTool) {
  return new Promise((resolve, reject) => {
    const worker = new Worker('src/workers/fasta-worker.js');
    const records = [];
    const t0 = performance.now();
    showError(`Parsing ${file.name}… 0 contigs so far`);

    worker.onmessage = (e) => {
      const msg = e.data;
      if (msg.type === 'contig') {
        records.push(msg.record);
      } else if (msg.type === 'progress') {
        showError(`Parsing ${file.name}… ${msg.contigsSoFar.toLocaleString()} contigs so far`);
      } else if (msg.type === 'done') {
        const elapsed = ((performance.now() - t0) / 1000).toFixed(1);
        showError(null);
        renderContigTable(records, binTablesByTool);
        document.getElementById('hMeta').textContent =
          `${msg.summary.contigCount.toLocaleString()} contigs · ${msg.summary.totalLength.toLocaleString()} bp · parsed in ${elapsed}s`;
        document.getElementById('hTitle').textContent = file.name;
        worker.terminate();
        resolve();
      } else if (msg.type === 'error') {
        showError(`Failed to parse ${file.name}: ${msg.message}`);
        worker.terminate();
        reject(new Error(msg.message));
      }
    };
    worker.onerror = (err) => {
      showError(`Failed to parse ${file.name}: ${err.message}`);
      worker.terminate();
      reject(err);
    };
    worker.postMessage({ file });
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
 * that looks like a contig->bin table (any number — Phase 5's cross-tool
 * reconciliation needs two or more, Phase 4's per-bin summary works with
 * just one) and parse it. The tool label per table comes from its
 * filename — content alone has no field naming which binning tool
 * produced it — de-duplicated if two files share a base name.
 * @returns {Promise<Map<string, {contigId:string,binId:string}[]>|null>}
 */
async function findAllBinTables(otherFiles) {
  const { sniff } = window.ClannMAG.sniff;
  const { parseContigBinTable } = window.ClannMAG.contigBinTable;
  const binTablesByTool = new Map();
  for (const file of otherFiles) {
    const text = await file.text();
    if (sniff(text).format !== 'contig-bin-table') continue;
    let label = fileBaseName(file.name);
    let suffix = 2;
    while (binTablesByTool.has(label)) label = `${fileBaseName(file.name)} (${suffix++})`;
    binTablesByTool.set(label, parseContigBinTable(text));
  }
  return binTablesByTool.size ? binTablesByTool : null;
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
    const binTablesByTool = otherFiles.length ? await findAllBinTables(otherFiles) : null;

    try {
      await loadAssembly(assemblyFile, binTablesByTool);
    } catch {
      // already surfaced via showError inside loadAssembly
    }
  });
}

document.addEventListener('DOMContentLoaded', () => {
  initTheme();
  initFilePicker();
});
})();
