(function () {
  'use strict';

// CSV serializers for Phase 9's export outputs (brief §Export): the
// revised contig->bin assignment table and the per-MAG/per-bin summary
// table. Pure string-building only — no Blob/download glue, which lives
// in app.js (browser-only, matches the project's pure-logic/DOM-glue
// split elsewhere).

function escapeCsvField(value) {
  const s = String(value);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

function toCsv(headerRow, dataRows) {
  const lines = [headerRow, ...dataRows].map((row) => row.map(escapeCsvField).join(','));
  return lines.join('\n') + '\n';
}

/** @param {{contigId:string, binId:string}[]} rows - working-assignment.js's assignmentToRows output */
function assignmentToCsv(rows) {
  return toCsv(['contig_id', 'bin_id'], rows.map((r) => [r.contigId, r.binId]));
}

/**
 * @param {object[]} summaries - bin-summary.js's computeBinSummaries output
 *   (binId, contigCount, totalLength, n50, meanGc, completeness, redundancy, mimagTier)
 * @param {Map<string, number>} [agreementFractionByBinId] - optional
 *   tool-agreement fraction per bin (e.g. mean contigAgreement.agreementFraction
 *   over a MAG's own contigs); omitted entries render blank
 */
function binSummaryToCsv(summaries, agreementFractionByBinId) {
  const header = ['bin_id', 'contig_count', 'total_length_bp', 'n50_bp', 'mean_gc_pct', 'completeness_pct', 'redundancy_pct', 'tool_agreement_fraction', 'mimag_tier'];
  const rows = summaries.map((s) => {
    const agreement = agreementFractionByBinId ? agreementFractionByBinId.get(s.binId) : undefined;
    return [
      s.binId, s.contigCount, s.totalLength, s.n50,
      s.meanGc != null ? (s.meanGc * 100).toFixed(2) : '',
      s.completeness.toFixed(2), s.redundancy.toFixed(2),
      agreement != null ? agreement.toFixed(3) : '',
      s.mimagTier,
    ];
  });
  return toCsv(header, rows);
}

const exportsObj = { assignmentToCsv, binSummaryToCsv };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.exportCsv = exportsObj;
}
})();
