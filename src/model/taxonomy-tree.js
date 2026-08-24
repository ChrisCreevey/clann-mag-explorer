(function () {
  'use strict';

// Flat, typed-array-shaped taxonomy tree keyed by taxid, ported from
// clann-edna-explorer/src/model/taxonomy-tree.js (see docs/phase1-investigation.md
// §2 for what changed vs. the original). Reused verbatim here: the class
// itself is generic over how nodes get populated.
//
// Two different population paths use this same tree in this app:
//  - build/03-taxonomy.js populates it offline from an NCBI taxdump
//    (nodes.dmp/names.dmp) to produce the shipped taxid->lineage table —
//    a flat parent-pointer walk, not indentation-based like a .breport.
//  - src/parsers/breport.js (optional per-contig Kraken2 input) populates
//    it in-browser from a .breport's indentation structure, unchanged from
//    the eDNA Explorer.
//
// Nodes are appended in first-seen order and addressed by a dense integer
// index; `taxidToIndex` maps taxid -> index for O(1) lookup.

const RANK_SUFFIX_RE = /^([A-Za-z])([0-9]*)$/;

function canonicalRank(rawRankCode) {
  const m = RANK_SUFFIX_RE.exec(rawRankCode.trim());
  if (!m) return { letter: rawRankCode.trim(), sub: 0 };
  return { letter: m[1], sub: m[2] ? Number(m[2]) : 0 };
}

class TaxonomyTree {
  constructor() {
    this.taxid = [];
    this.name = [];
    this.rankLetter = [];
    this.rankSub = [];
    this.depth = [];
    this.parentIndex = [];
    this.perSample = []; // perSample[index] = Map<key, {...}> — key is a contigId or sampleId depending on caller
    this.taxidToIndex = new Map();
  }

  get size() {
    return this.taxid.length;
  }

  getOrCreateNode(taxid, name, rawRankCode, depth, parentTaxid) {
    let idx = this.taxidToIndex.get(taxid);
    if (idx !== undefined) return idx;

    const { letter, sub } = canonicalRank(rawRankCode);
    idx = this.taxid.length;
    this.taxid.push(taxid);
    this.name.push(name);
    this.rankLetter.push(letter);
    this.rankSub.push(sub);
    this.depth.push(depth);
    this.parentIndex.push(
      parentTaxid === null ? -1 : this.taxidToIndex.get(parentTaxid) ?? -1
    );
    this.perSample.push(new Map());
    this.taxidToIndex.set(taxid, idx);
    return idx;
  }

  setSampleCounts(taxid, key, counts) {
    const idx = this.taxidToIndex.get(taxid);
    if (idx === undefined) {
      throw new Error(`setSampleCounts: unknown taxid ${taxid}`);
    }
    const existing = this.perSample[idx].get(key) || {};
    this.perSample[idx].set(key, { ...existing, ...counts });
  }

  getSampleCounts(taxid, key) {
    const idx = this.taxidToIndex.get(taxid);
    if (idx === undefined) return undefined;
    return this.perSample[idx].get(key);
  }

  node(taxid) {
    const idx = this.taxidToIndex.get(taxid);
    if (idx === undefined) return null;
    return {
      taxid: this.taxid[idx],
      name: this.name[idx],
      rank: this.rankSub[idx] ? `${this.rankLetter[idx]}${this.rankSub[idx]}` : this.rankLetter[idx],
      depth: this.depth[idx],
      parentTaxid: this.parentIndex[idx] === -1 ? null : this.taxid[this.parentIndex[idx]],
    };
  }

  isLeafRank(taxid, leafRankLetter = 'S') {
    const idx = this.taxidToIndex.get(taxid);
    if (idx === undefined) return false;
    return this.rankLetter[idx] === leafRankLetter && this.rankSub[idx] === 0;
  }

  // --- New vs. eDNA Explorer: no LCA function existed there to port (see
  // docs/phase1-investigation.md §2). Needed here for the marker-gene
  // taxonomic consistency check (brief §Marker-gene identification module).

  /** Root-to-node path of indices, root first. */
  pathIndices(taxid) {
    const idx = this.taxidToIndex.get(taxid);
    if (idx === undefined) return [];
    const path = [];
    let cur = idx;
    while (cur !== -1) {
      path.push(cur);
      cur = this.parentIndex[cur];
    }
    return path.reverse();
  }

  /** Lowest common ancestor taxid of a list of taxids, or null if none share a root. */
  lca(taxids) {
    const paths = taxids
      .map((t) => this.pathIndices(t))
      .filter((p) => p.length > 0);
    if (paths.length === 0) return null;
    let i = 0;
    const minLen = Math.min(...paths.map((p) => p.length));
    while (i < minLen && paths.every((p) => p[i] === paths[0][i])) i++;
    if (i === 0) return null;
    return this.taxid[paths[0][i - 1]];
  }
}

const taxonomyTreeExports = { TaxonomyTree, canonicalRank };
if (typeof module !== 'undefined' && module.exports) module.exports = taxonomyTreeExports;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.taxonomyTree = taxonomyTreeExports;
}
})();
