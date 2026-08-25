(function () {
  'use strict';

// Parses a coverage depth table: contig ID plus one depth column per
// sample (brief §Inputs). Two real shapes handled:
//  - MetaBAT2's `jgi_summarize_bam_contig_depths` output: a header row
//    `contigName  contigLen  totalAvgDepth  sample1.bam  sample1.bam-var
//    sample2.bam  sample2.bam-var ...` — depth and variance columns
//    interleaved per sample. Only the depth columns become samples; the
//    "-var" columns and the derived contigLen/totalAvgDepth columns are
//    dropped (variance isn't a per-sample depth, and length is already
//    computed from the assembly itself).
//  - A generic table: contig ID + one numeric column per sample, with or
//    without a header row.

const KNOWN_NON_SAMPLE_COLUMNS = new Set(['contigname', 'contig', 'contig_id', 'contiglen', 'totalavgdepth']);

function detectDelimiter(line) {
  const tabCount = (line.match(/\t/g) || []).length;
  const commaCount = (line.match(/,/g) || []).length;
  return tabCount >= commaCount ? '\t' : ',';
}

function looksLikeHeaderRow(fields) {
  // A header row's non-first cells won't parse as numbers; a data row's will.
  return fields.slice(1).some((f) => f.trim() !== '' && Number.isNaN(Number(f)));
}

/**
 * @param {string} text
 * @returns {{sampleNames: string[], rows: {contigId: string, depths: number[]}[]}}
 */
function parseCoverageTable(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
  if (lines.length === 0) return { sampleNames: [], rows: [] };

  const delim = detectDelimiter(lines[0]);
  const firstFields = lines[0].split(delim);
  const hasHeader = looksLikeHeaderRow(firstFields);

  let sampleNames;
  let sampleColumnIndices;
  if (hasHeader) {
    sampleColumnIndices = [];
    sampleNames = [];
    for (let i = 1; i < firstFields.length; i++) {
      const label = firstFields[i].trim();
      if (KNOWN_NON_SAMPLE_COLUMNS.has(label.toLowerCase()) || /-var$/i.test(label)) continue;
      sampleColumnIndices.push(i);
      sampleNames.push(label);
    }
  } else {
    sampleColumnIndices = firstFields.slice(1).map((_, i) => i + 1);
    sampleNames = sampleColumnIndices.map((_, i) => `sample${i + 1}`);
  }

  const dataLines = hasHeader ? lines.slice(1) : lines;
  const rows = [];
  for (const line of dataLines) {
    const fields = line.split(delim);
    const contigId = fields[0].trim();
    if (!contigId) continue;
    const depths = sampleColumnIndices.map((i) => Number(fields[i]));
    rows.push({ contigId, depths });
  }

  return { sampleNames, rows };
}

/** Cheap content-based check for sniff.js. */
function looksLikeCoverageTable(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0).slice(0, 10);
  if (lines.length === 0) return false;
  const delim = detectDelimiter(lines[0]);
  const firstFields = lines[0].split(delim);
  if (firstFields.length < 2) return false;
  const dataLines = looksLikeHeaderRow(firstFields) ? lines.slice(1) : lines;
  if (dataLines.length === 0) return false;
  return dataLines.every((line) => {
    const fields = line.split(delim);
    if (fields.length < 2 || !fields[0].trim()) return false;
    return fields.slice(1).every((f) => f.trim() === '' || !Number.isNaN(Number(f)));
  });
}

const exportsObj = { parseCoverageTable, looksLikeCoverageTable };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.coverageTable = exportsObj;
}
})();
