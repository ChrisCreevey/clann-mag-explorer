(function () {
  'use strict';

// App shell. Loads an assembly FASTA (streamed through fasta-index.js in a
// Worker, with marker-gene search) and, per Phase 4, an optional single
// contig->bin table alongside it — content-sniffed from whatever else was
// selected, not by filename. Cross-tool reconciliation (multiple bin
// tables at once), filters, and reassignment land in later phases.

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

function renderBinSummaryCard(records, binAssignments) {
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
      <h3>Bin summaries</h3>
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

function renderContigTable(records, binAssignments) {
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
    ${binAssignments ? renderBinSummaryCard(records, binAssignments) : ''}
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
function loadAssembly(file, binAssignments) {
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
        renderContigTable(records, binAssignments);
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

/**
 * Among the files NOT identified as the assembly, find the first one that
 * content-sniffs as a contig->bin table and parse it. Multiple bin tables
 * (cross-tool reconciliation) is Phase 5 scope — for now, first match
 * wins and the rest are ignored silently, matching how the assembly pick
 * itself already works.
 */
async function findBinAssignments(otherFiles) {
  const { sniff } = window.ClannMAG.sniff;
  const { parseContigBinTable } = window.ClannMAG.contigBinTable;
  for (const file of otherFiles) {
    const text = await file.text();
    if (sniff(text).format === 'contig-bin-table') {
      return parseContigBinTable(text);
    }
  }
  return null;
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
    const binAssignments = otherFiles.length ? await findBinAssignments(otherFiles) : null;

    try {
      await loadAssembly(assemblyFile, binAssignments);
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
