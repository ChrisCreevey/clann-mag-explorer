(function () {
  'use strict';

// Content-based format detection, matching clann-edna-explorer's convention:
// files are identified by their content, not filename or extension.
// FASTA/gzip detection lives separately in app.js (byte-level peek, needed
// before any text decoding happens) — this covers the plain-text tabular
// inputs.

const { parseBreportLine } = (typeof module !== 'undefined' && module.exports)
  ? require('./breport')
  : self.ClannMAG.breport;
const { looksLikeContigBinTable } = (typeof module !== 'undefined' && module.exports)
  ? require('./contig-bin-table')
  : self.ClannMAG.contigBinTable;
const { looksLikeCoverageTable } = (typeof module !== 'undefined' && module.exports)
  ? require('./coverage-table')
  : self.ClannMAG.coverageTable;
const { looksLikeKraken2ContigOutput } = (typeof module !== 'undefined' && module.exports)
  ? require('./kraken2-contigs')
  : self.ClannMAG.kraken2Contigs;

function looksLikeBreport(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0).slice(0, 20);
  if (lines.length === 0) return false;
  return lines.every((line) => parseBreportLine(line) !== null);
}

function columnCount(text) {
  const firstLine = text.split(/\r?\n/).find((l) => l.trim().length > 0) || '';
  const tabCount = (firstLine.match(/\t/g) || []).length;
  const commaCount = (firstLine.match(/,/g) || []).length;
  return 1 + Math.max(tabCount, commaCount);
}

/**
 * @returns {{format: 'breport'|'kraken2-contigs'|'contig-bin-table'|'coverage-table'|'unknown', reason?: string}}
 */
function sniff(text) {
  // Order matters where shapes could otherwise collide:
  //  - breport (exactly 6 tab fields, strict per-field shape) and
  //    kraken2-contigs (C/U flag in field 1) are both distinctive enough
  //    to check first with no ambiguity risk.
  //  - A plain 2-column table is ambiguous between contig-bin-table
  //    (contig, binId) and a single-sample coverage-table (contig, depth)
  //    when the second column happens to be numeric either way. Per the
  //    documented heuristic in docs/phase1-investigation.md (brief's own
  //    "two-column" vs "per-sample depth columns" framing), 2 columns
  //    defaults to contig-bin-table; 3+ columns (multi-sample depth, or
  //    MetaBAT2's contigName/contigLen/totalAvgDepth/... header) is what
  //    distinguishes a real coverage table.
  if (looksLikeBreport(text)) return { format: 'breport' };
  if (looksLikeKraken2ContigOutput(text)) return { format: 'kraken2-contigs' };
  if (columnCount(text) >= 3 && looksLikeCoverageTable(text)) return { format: 'coverage-table' };
  if (looksLikeContigBinTable(text)) return { format: 'contig-bin-table' };
  if (looksLikeCoverageTable(text)) return { format: 'coverage-table' };
  return { format: 'unknown', reason: 'did not match any recognised tabular input shape' };
}

const sniffExports = { sniff };
if (typeof module !== 'undefined' && module.exports) module.exports = sniffExports;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.sniff = sniffExports;
}
})();
