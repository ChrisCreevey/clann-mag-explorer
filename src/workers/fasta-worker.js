'use strict';

// Runs the streaming FASTA parse off the main thread (Phase 2 performance
// follow-up — see docs/phase1-investigation.md "Phase 2 findings"). This
// doesn't reduce the total work, it just stops a multi-second-to-minutes
// parse from freezing the tab: the main thread posts a File, this worker
// streams it and posts per-contig records + periodic progress back, so the
// UI can render a live progress bar and stay responsive throughout.
//
// Phase 3: also runs the marker-gene search (src/model/marker-genes.js)
// against each contig's six-frame translation, reusing the same
// translation contig-stats.js already computed (record.frames) rather
// than recomputing it. Marker-gene assets (~27MB) start loading as soon
// as this worker starts, in parallel with whatever the user does before
// picking a file, so they're often already resolved by the time parsing
// begins; if they're not (or fail to load — network hiccup, assets
// missing), marker search is skipped rather than blocking or failing the
// whole load, per the brief's framing of this as an optional module.
//
// importScripts (not <script> tags) since Workers don't share the page's
// DOM/script context; relative to this file's own location, not the page's.
importScripts(
  '../model/dna-codes.js', '../model/translate.js', '../model/contig-stats.js', '../model/fasta-index.js',
  '../model/reduced-alphabet.js', '../model/blosum62.js', '../model/marker-genes.js'
);

const { streamFasta } = self.ClannMAG.fastaIndex;
const { loadMarkerGeneAssets, searchContigForMarkers } = self.ClannMAG.markerGenes;

const PROGRESS_EVERY = 200; // contigs between progress pings, not every single one

// Resolve relative to this worker script's own URL, not the page's, so it
// works regardless of the site's deployment path (e.g. GitHub Pages
// project subpaths).
const DATA_URL = new URL('../../data/', self.location).href;

const markerAssetsPromise = loadMarkerGeneAssets(DATA_URL).catch((err) => {
  self.postMessage({ type: 'marker-assets-error', message: err.message });
  return null;
});

self.onmessage = async (e) => {
  const { file } = e.data;
  const markerAssets = await markerAssetsPromise;
  let count = 0;
  try {
    const summary = await streamFasta(file, (record) => {
      count++;
      if (markerAssets && record.frames) {
        record.markerHits = searchContigForMarkers(record.frames, markerAssets);
      }
      delete record.frames; // transient — only needed for the marker search just above
      self.postMessage({ type: 'contig', record });
      if (count % PROGRESS_EVERY === 0) {
        self.postMessage({ type: 'progress', contigsSoFar: count });
      }
    });
    self.postMessage({ type: 'done', summary });
  } catch (err) {
    self.postMessage({ type: 'error', message: err.message });
  }
};
