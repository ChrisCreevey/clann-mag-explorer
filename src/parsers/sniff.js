(function () {
  'use strict';

// Content-based format detection, matching clann-edna-explorer's convention:
// files are identified by their content, not filename or extension.
// FASTA/gzip detection lives separately in app.js (byte-level peek, needed
// before any text decoding happens) — this covers the plain-text tabular
// inputs. Only distinguishes formats this repo can actually parse yet:
// coverage-table.js is still a stub, so sniff() doesn't route to it until
// that lands, to avoid detecting a format it can't handle.

const { parseBreportLine } = (typeof module !== 'undefined' && module.exports)
  ? require('./breport')
  : self.ClannMAG.breport;
const { looksLikeContigBinTable } = (typeof module !== 'undefined' && module.exports)
  ? require('./contig-bin-table')
  : self.ClannMAG.contigBinTable;

function looksLikeBreport(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0).slice(0, 20);
  if (lines.length === 0) return false;
  return lines.every((line) => parseBreportLine(line) !== null);
}

/** @returns {{format: 'breport'|'contig-bin-table'|'unknown', reason?: string}} */
function sniff(text) {
  if (looksLikeBreport(text)) return { format: 'breport' };
  if (looksLikeContigBinTable(text)) return { format: 'contig-bin-table' };
  return { format: 'unknown', reason: 'did not match any recognised tabular input shape' };
}

const sniffExports = { sniff };
if (typeof module !== 'undefined' && module.exports) module.exports = sniffExports;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.sniff = sniffExports;
}
})();
