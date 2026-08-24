(function () {
  'use strict';

// Marker-gene identification module (brief §Marker-gene identification
// module): loads the offline-built seed index, reference sequences, and
// taxid->lineage table (produced by build/01-cluster.js, build/02-index.js,
// build/03-taxonomy.js — see docs/phase1-investigation.md §4-5), and runs
// the seed-and-extend search against each contig's six-frame translation.
// Stubbed pending Phase 3, which itself is blocked on the build/ pipeline
// existing and its assets being calibrated.

async function loadMarkerGeneAssets(dataUrl) {
  throw new Error('loadMarkerGeneAssets: not implemented yet — build/ pipeline outputs required first');
}

function searchContigForMarkers(sixFrameTranslations, assets) {
  throw new Error('searchContigForMarkers: not implemented yet');
}

const exportsObj = { loadMarkerGeneAssets, searchContigForMarkers };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof window !== 'undefined') {
  window.ClannMAG = window.ClannMAG || {};
  window.ClannMAG.markerGenes = exportsObj;
}
})();
