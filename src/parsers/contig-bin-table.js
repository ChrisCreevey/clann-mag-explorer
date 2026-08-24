(function () {
  'use strict';

// Parses a two-column tab/CSV contig->bin assignment table (one per binning
// tool). Stubbed pending Phase 1 investigation into real MetaBAT2/CONCOCT/
// MaxBin2 output shapes.

function parseContigBinTable(text) {
  throw new Error('parseContigBinTable: not implemented yet');
}

const exportsObj = { parseContigBinTable };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof window !== 'undefined') {
  window.ClannMAG = window.ClannMAG || {};
  window.ClannMAG.contigBinTable = exportsObj;
}
})();
