(function () {
  'use strict';

// Per-contig statistics computed on the streaming FASTA pass: length, GC%,
// GC skew, tetranucleotide composition signature, coding-density estimate.
// Stubbed pending Phase 2.

function computeContigStats(sequence) {
  throw new Error('computeContigStats: not implemented yet');
}

const exportsObj = { computeContigStats };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof window !== 'undefined') {
  window.ClannMAG = window.ClannMAG || {};
  window.ClannMAG.contigStats = exportsObj;
}
})();
