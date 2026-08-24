(function () {
  'use strict';

// App shell. Phase 2 wires up just enough to load a single assembly FASTA
// end to end: pick a file, detect it looks like FASTA (gzip-aware peek,
// content-based per the suite's convention — a real sniff.js covering
// contig->bin/coverage tables too is still Phase 1/4 work), stream it
// through fasta-index.js, and render the per-contig stats table. Bin
// loading, filters, and reassignment land in later phases.

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

function renderContigTable(records) {
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
function loadAssembly(file) {
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
        renderContigTable(records);
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

function initFilePicker() {
  const input = document.getElementById('folder-input');
  const openButtons = [document.getElementById('uploadBtn'), document.getElementById('emptyOpen')];
  openButtons.forEach((btn) => btn && btn.addEventListener('click', () => input.click()));
  input.addEventListener('change', async () => {
    const files = [...input.files];
    if (files.length === 0) return;

    // Contig->bin / coverage / breport loading is Phase 4+; for now, load
    // the first file that content-sniffs as a FASTA and ignore the rest,
    // rather than silently doing nothing.
    let assemblyFile = null;
    for (const file of files) {
      if (await looksLikeFasta(file)) { assemblyFile = file; break; }
    }
    if (!assemblyFile) {
      showError('No file among your selection looks like a FASTA assembly (checked for a leading ">" or gzip magic bytes). Bin/coverage table loading isn’t implemented yet.');
      return;
    }
    try {
      await loadAssembly(assemblyFile);
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
