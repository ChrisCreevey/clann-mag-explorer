(function () {
  'use strict';

// MAG-level filtering, separate from contig-level filtering (filters.js):
// the user asked for two distinct filter axes — one over per-contig
// properties, one over the cross-tool reconciliation table's own columns
// (Putative MAG, Core, Disputed, Completeness, Redundancy, Tier, and
// which tool(s) support it). Selecting MAGs here is meant to reduce every
// MAG-scoped table and graph (the reconciliation table, the agreement
// network, the QC section) down to just the chosen MAGs — a different
// question from "which contigs pass a property filter", so it gets its
// own filter state and its own left-pane section rather than being folded
// into filters.js.

function defaultMagFilters() {
  return {
    magIdSearch: '',
    tiers: { high: true, medium: true, low: true },
    coreMin: null, coreMax: null,
    disputedMin: null, disputedMax: null,
    completenessMin: null, completenessMax: null,
    redundancyMin: null, redundancyMax: null,
    supportedByTool: '', // '' = any tool
  };
}

function passesRange(value, min, max) {
  if (min != null && value < min) return false;
  if (max != null && value > max) return false;
  return true;
}

/**
 * @param {{magId:string, coreCount:number, disputedCount:number,
 *   completeness:number, redundancy:number, tier:'high'|'medium'|'low',
 *   tools:string[]}[]} magSummaryData - one row per putative MAG, the same
 *   numbers the reconciliation table itself renders
 * @param {object} filters - see defaultMagFilters()
 * @returns {object[]} the subset of magSummaryData passing every filter
 */
function applyMagFilters(magSummaryData, filters) {
  const needle = (filters.magIdSearch || '').trim().toLowerCase();
  return magSummaryData.filter((m) => {
    if (needle && !m.magId.toLowerCase().includes(needle)) return false;
    if (filters.tiers && filters.tiers[m.tier] === false) return false;
    if (!passesRange(m.coreCount, filters.coreMin, filters.coreMax)) return false;
    if (!passesRange(m.disputedCount, filters.disputedMin, filters.disputedMax)) return false;
    if (!passesRange(m.completeness, filters.completenessMin, filters.completenessMax)) return false;
    if (!passesRange(m.redundancy, filters.redundancyMin, filters.redundancyMax)) return false;
    if (filters.supportedByTool && !m.tools.includes(filters.supportedByTool)) return false;
    return true;
  });
}

const exportsObj = { defaultMagFilters, applyMagFilters };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.magFilters = exportsObj;
}
})();
