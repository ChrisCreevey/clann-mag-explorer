(function () {
  'use strict';

// Parser for Kraken-style .breport files, ported verbatim from
// clann-edna-explorer/src/parsers/breport.js — one row per taxonomy tree
// node, hierarchy encoded via leading-space indentation on the name column
// (2 spaces per depth level). No header row.
//
// Columns: pct_reads_rooted, reads_covered(clade), reads_assigned(direct),
// rank_code, taxid, indented_name.
//
// Here the report is run against the assembly's contigs rather than reads
// (see brief §Inputs), so callers pass a contig ID in place of eDNA
// Explorer's sampleId — TaxonomyTree.setSampleCounts is generic over that
// key, so no change to the tree/parser API is needed for that swap.
//
// Streams line-by-line rather than building an intermediate array of every
// line.

function parseBreportLine(line) {
  const cols = line.split('\t');
  if (cols.length !== 6) return null;
  const [pctStr, cladeStr, directStr, rankCode, taxidStr, rawName] = cols;

  const leadingSpaces = rawName.length - rawName.replace(/^ +/, '').length;
  const depth = leadingSpaces / 2;

  return {
    pctOfTotal: Number(pctStr),
    cladeReads: Number(cladeStr),
    directReads: Number(directStr),
    rankCode: rankCode.trim(),
    taxid: Number(taxidStr),
    name: rawName.trim(),
    depth,
  };
}

/**
 * Parse .breport text and attach nodes/counts into the given shared
 * TaxonomyTree for `contigId`. Returns summary stats for the contig.
 *
 * @param {string} text - full file contents
 * @param {import('../model/taxonomy-tree').TaxonomyTree} tree
 * @param {string} contigId
 */
function parseBreport(text, tree, contigId) {
  const lines = text.split(/\r?\n/);
  // stack[d] = taxid of the current open ancestor at depth d
  const stack = [];
  let rowCount = 0;
  let rootCladeReads = 0;

  for (const rawLine of lines) {
    if (rawLine.trim() === '') continue;
    const row = parseBreportLine(rawLine);
    if (!row) continue;
    rowCount++;

    const parentTaxid = row.depth === 0 ? null : stack[row.depth - 1] ?? null;
    tree.getOrCreateNode(row.taxid, row.name, row.rankCode, row.depth, parentTaxid);
    tree.setSampleCounts(row.taxid, contigId, {
      cladeReads: row.cladeReads,
      directReads: row.directReads,
      pctOfTotal: row.pctOfTotal,
    });

    stack[row.depth] = row.taxid;
    stack.length = row.depth + 1;

    if (row.depth === 0) rootCladeReads = row.cladeReads;
  }

  return {
    rowCount,
    totalReads: rootCladeReads,
  };
}

const breportExports = { parseBreport, parseBreportLine };
if (typeof module !== 'undefined' && module.exports) module.exports = breportExports;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.breport = breportExports;
}
})();
