(function () {
  'use strict';

// Parses a two-column tab/CSV contig->bin assignment table (brief
// §Inputs), one per binning tool. Matches DAS_Tool's widely-used
// "scaffolds2bin.tsv" convention (its own Fasta_to_Contig2Bin.sh helper
// produces exactly this: contig<TAB>bin, no header) and CONCOCT's
// clustering_gt1000.csv (contig_id,cluster_id, comma-delimited, WITH a
// header row) — so both delimiter and header presence are detected from
// content, not assumed, per the brief's content-based-detection
// convention used throughout this suite.

const HEADER_FIRST_COLUMN_NAMES = new Set([
  'contig', 'contigid', 'contig_id', 'contigname', 'contig_name',
  'scaffold', 'scaffoldid', 'scaffold_id', 'sequence', 'sequence_id', 'sequenceid', 'name',
]);

function detectDelimiter(line) {
  const tabCount = (line.match(/\t/g) || []).length;
  const commaCount = (line.match(/,/g) || []).length;
  return tabCount >= commaCount ? '\t' : ',';
}

function looksLikeHeaderRow(fields) {
  return HEADER_FIRST_COLUMN_NAMES.has(fields[0].trim().toLowerCase());
}

/**
 * @param {string} text - full file contents
 * @returns {{contigId: string, binId: string}[]}
 */
function parseContigBinTable(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
  if (lines.length === 0) return [];

  const delim = detectDelimiter(lines[0]);
  const firstFields = lines[0].split(delim);
  const startIdx = looksLikeHeaderRow(firstFields) ? 1 : 0;

  const assignments = [];
  for (let i = startIdx; i < lines.length; i++) {
    const fields = lines[i].split(delim);
    if (fields.length < 2) continue;
    const contigId = fields[0].trim();
    const binId = fields[1].trim();
    if (!contigId || !binId) continue;
    assignments.push({ contigId, binId });
  }
  return assignments;
}

/**
 * Cheap content-based check for sniff.js: does this look like a two-column
 * contig->bin table at all? Not a full parse — just "would parsing this
 * produce at least one plausible row without the shape falling apart."
 */
function looksLikeContigBinTable(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0).slice(0, 20);
  if (lines.length === 0) return false;
  const delim = detectDelimiter(lines[0]);
  const dataLines = looksLikeHeaderRow(lines[0].split(delim)) ? lines.slice(1) : lines;
  if (dataLines.length === 0) return false;
  return dataLines.every((line) => {
    const fields = line.split(delim);
    return fields.length === 2 && fields[0].trim().length > 0 && fields[1].trim().length > 0;
  });
}

const exportsObj = { parseContigBinTable, looksLikeContigBinTable };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.contigBinTable = exportsObj;
}
})();
