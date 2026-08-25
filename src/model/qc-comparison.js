(function () {
  'use strict';

// Pure helpers for Phase 8's "comparison and QC across the full set"
// (brief): picking which bins to show in the good-vs-bad side-by-side
// view, and shaping putative-MAG data for the completeness/contamination
// scatter. Kept separate from bin-summary.js/bin-reconciliation.js since
// those compute the underlying numbers — this module only picks and
// reshapes for display, no new statistics.

/**
 * Picks the best- and worst-quality bin from a set of bin summaries
 * (bin-summary.js's computeBinSummaries output) for direct visual
 * contrast — brief's "tight, well-separated cluster vs.
 * scattered/fragmented one". Score = completeness minus redundancy, so a
 * bin that is both complete and low-contamination ranks best; only bins
 * with at least 2 contigs are eligible (a single-contig bin has no
 * cluster shape to contrast against).
 * @param {object[]} summaries - bin-summary.js summaries (each has
 *   `binId`, `contigCount`, `completeness`, `redundancy`)
 * @returns {{good: object|null, bad: object|null}} the chosen summary
 *   objects (or null if fewer than 2 eligible bins exist)
 */
function pickComparisonBins(summaries) {
  const eligible = summaries.filter((s) => s.contigCount >= 2);
  if (eligible.length < 2) return { good: null, bad: null };

  const scored = eligible.map((s) => ({ summary: s, score: s.completeness - s.redundancy }));
  scored.sort((a, b) => b.score - a.score);
  return { good: scored[0].summary, bad: scored[scored.length - 1].summary };
}

/**
 * Reshapes Phase 5's putative MAGs plus their bin-summary stats into
 * scatter-ready points for the completeness-vs-contamination view (brief:
 * "coloured by which tool(s) support each one").
 * @param {object[]} magSummaries - bin-summary.js summaries computed over
 *   an assignment keyed by magId (one point per putative MAG)
 * @param {{magId:string, members:{tool:string}[]}[]} putativeMags
 * @returns {{id:string, x:number, y:number, colorKey:string}[]} x =
 *   completeness, y = redundancy (contamination proxy), colorKey = the
 *   sorted, comma-joined list of supporting tool names
 */
function buildMagQcPoints(magSummaries, putativeMags) {
  const toolsByMagId = new Map(
    putativeMags.map((m) => [m.magId, [...new Set(m.members.map((mem) => mem.tool))].sort().join(', ')])
  );
  return magSummaries.map((s) => ({
    id: s.binId, x: s.completeness, y: s.redundancy, colorKey: toolsByMagId.get(s.binId) || 'unknown',
  }));
}

const exportsObj = { pickComparisonBins, buildMagQcPoints };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.qcComparison = exportsObj;
}
})();
