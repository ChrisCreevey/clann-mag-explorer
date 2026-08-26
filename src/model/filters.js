(function () {
  'use strict';

// Left-pane filtering (brief §Left pane: filters and search): length, GC%,
// coding density range filters, bin (from any loaded tool) or "unbinned",
// cross-tool agreement fraction, and a contig/bin ID search — all combined
// with AND semantics. Pure data-in/data-out so it's usable identically from
// app.js (against live records) and from tests, matching outliers.js/
// bin-reconciliation.js's separation of pure model from DOM wiring.

const BIN_FILTER_SEP = '␟'; // unlikely to collide with a real tool/bin name

function defaultFilters() {
  return {
    lengthMin: null, lengthMax: null,
    gcMin: null, gcMax: null,
    codingDensityMin: null, codingDensityMax: null,
    binFilter: '', // '' = any, '__unbinned__', or `${tool}${BIN_FILTER_SEP}${binId}`
    maxAgreementPercent: null, // null = no agreement filter
    searchText: '',
  };
}

/** @param {Map<string, {contigId:string, binId:string}[]>|null} binTablesByTool */
function buildBinIndex(binTablesByTool) {
  const index = new Map(); // tool -> Map<contigId, binId>
  if (!binTablesByTool) return index;
  for (const [tool, assignments] of binTablesByTool) {
    const m = new Map();
    for (const { contigId, binId } of assignments) m.set(contigId, binId);
    index.set(tool, m);
  }
  return index;
}

function listBinFilterOptions(binIndex) {
  const options = [];
  for (const [tool, byContig] of binIndex) {
    const binIds = [...new Set(byContig.values())].sort((a, b) => a.localeCompare(b));
    for (const binId of binIds) options.push({ value: `${tool}${BIN_FILTER_SEP}${binId}`, label: `${tool}: ${binId}`, tool, binId });
  }
  return options;
}

function passesRange(value, min, max) {
  if (min != null && value < min) return false;
  if (max != null && value > max) return false;
  return true;
}

function matchesBinFilter(contigId, binFilter, binIndex) {
  if (!binFilter) return true;
  if (binFilter === '__unbinned__') {
    for (const byContig of binIndex.values()) if (byContig.has(contigId)) return false;
    return true;
  }
  const sep = binFilter.indexOf(BIN_FILTER_SEP);
  if (sep < 0) return true;
  const tool = binFilter.slice(0, sep);
  const binId = binFilter.slice(sep + 1);
  const byContig = binIndex.get(tool);
  return !!byContig && byContig.get(contigId) === binId;
}

function matchesSearch(record, searchText, binIndex) {
  const needle = (searchText || '').trim().toLowerCase();
  if (needle === '') return true;
  if (record.id.toLowerCase().includes(needle)) return true;
  for (const byContig of binIndex.values()) {
    const binId = byContig.get(record.id);
    if (binId && binId.toLowerCase().includes(needle)) return true;
  }
  return false;
}

/**
 * @param {object[]} records
 * @param {object} filters - see defaultFilters()
 * @param {{binIndex?: Map, agreementByContigId?: Map<string,number>}} ctx
 */
function applyFilters(records, filters, ctx = {}) {
  const binIndex = ctx.binIndex || new Map();
  const agreementByContigId = ctx.agreementByContigId || new Map();
  return records.filter((r) => {
    if (!passesRange(r.length, filters.lengthMin, filters.lengthMax)) return false;
    if (!passesRange(r.gcContent * 100, filters.gcMin, filters.gcMax)) return false;
    if (!passesRange(r.codingDensity * 100, filters.codingDensityMin, filters.codingDensityMax)) return false;
    if (!matchesBinFilter(r.id, filters.binFilter, binIndex)) return false;
    if (filters.maxAgreementPercent != null) {
      const agreement = agreementByContigId.get(r.id);
      if (agreement == null || agreement * 100 > filters.maxAgreementPercent) return false;
    }
    if (!matchesSearch(r, filters.searchText, binIndex)) return false;
    return true;
  });
}

const exportsObj = { defaultFilters, buildBinIndex, listBinFilterOptions, applyFilters, BIN_FILTER_SEP };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.filters = exportsObj;
}
})();
