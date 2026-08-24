(function () {
  'use strict';

// Content-based format detection, matching clann-edna-explorer's convention:
// files are identified by their content, not filename or extension.
// Stubbed pending Phase 1 calibration of contig->bin / coverage table shapes
// against real binning-tool output.

function sniff(text) {
  return { format: 'unknown', reason: 'sniffing not implemented yet' };
}

const sniffExports = { sniff };
if (typeof module !== 'undefined' && module.exports) module.exports = sniffExports;
if (typeof window !== 'undefined') {
  window.ClannMAG = window.ClannMAG || {};
  window.ClannMAG.sniff = sniffExports;
}
})();
