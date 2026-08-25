(function () {
  'use strict';

// Marker-gene taxonomic consistency (brief §Marker-gene identification
// module, "using provenance the search already produces for free"): each
// called marker gene carries a provenance taxID (its top-scoring
// reference's nearest known relative), so collecting those across a
// bin's markers and computing their lowest common ancestor gives a
// chimerism signal — a bin whose markers only agree at an unexpectedly
// coarse rank is a chimerism candidate; per-contig, whichever contig's
// markers sit furthest from the bin's consensus becomes an outlier
// signal, alongside composition/coverage distance (outliers.js) and
// marker redundancy (bin-summary.js's computeMarkerContributions).
//
// Built on taxonomy-tree.js's generic TaxonomyTree/lca() — see
// docs/phase1-investigation.md §2 for why that needed a new LCA function
// (the eDNA Explorer's tree had the data shape but no LCA) and §4 for why
// this module's lineage table (build/03-taxonomy.js's output, resolved
// against a real NCBI taxdump) is a different ingestion path than
// breport.js's indentation-derived tree.

const { TaxonomyTree } = (typeof module !== 'undefined' && module.exports)
  ? require('./taxonomy-tree')
  : self.ClannMAG.taxonomyTree;

/**
 * @param {{taxid:number, parentTaxid:number, rank:string, name:string}[]} lineageRows
 *   build/03-taxonomy.js's output (data/scg40-lineage.json)
 */
function buildTaxonomyTreeFromLineage(lineageRows) {
  const tree = new TaxonomyTree();
  const byTaxid = new Map(lineageRows.map((n) => [n.taxid, n]));
  const inserted = new Set();

  function insert(taxid) {
    if (inserted.has(taxid)) return;
    const row = byTaxid.get(taxid);
    if (!row) return; // shouldn't happen — build script includes every ancestor — but don't crash if it does
    const isRoot = row.parentTaxid === taxid;
    if (!isRoot) insert(row.parentTaxid);
    const parentTaxid = isRoot ? null : row.parentTaxid;
    const depth = isRoot ? 0 : tree.node(parentTaxid).depth + 1;
    tree.getOrCreateNode(taxid, row.name, row.rank, depth, parentTaxid);
    inserted.add(taxid);
  }

  for (const row of lineageRows) insert(row.taxid);
  return tree;
}

/** Fetch + build in one step, for callers that just want a ready tree. */
async function loadTaxonomyTree(dataUrl) {
  const lineageRows = await fetch(dataUrl + 'scg40-lineage.json').then((r) => r.json());
  return buildTaxonomyTreeFromLineage(lineageRows);
}

/**
 * @param {object[]} binContigs - per-contig records with `markerHits:
 *   {provenanceTaxId:number}[]`
 * @param {TaxonomyTree} tree
 * @returns {{consensusTaxId:number|null, consensusRank:string|null,
 *   consensusName:string|null, perContigDistance:Map<string, number|null>}}
 *   perContigDistance: a leave-one-out measure — how many ranks "up" the
 *   *rest of the bin's own* consensus gets pulled by adding this contig's
 *   markers back in. 0 means this contig's provenance is already
 *   consistent with (a descendant of) what the rest of the bin agrees
 *   on; higher means this contig is dragging the group toward a coarser,
 *   less specific common ancestor — a chimerism candidate. Note this is
 *   deliberately NOT "distance from the whole-bin consensus": the
 *   whole-bin LCA is by definition an ancestor of every individual
 *   contig's own LCA, so combining any one contig's taxIDs back into it
 *   is a mathematical no-op (always returns the same consensus) — caught
 *   before shipping by checking the math, not by a failing test. Leaving
 *   this contig OUT first, so the comparison point isn't already forced
 *   to already include it, is what makes the signal real.
 */
function computeBinTaxonomicConsistency(binContigs, tree) {
  const perContigOwnTaxIds = binContigs.map(
    (contig) => (contig.markerHits || []).map((h) => h.provenanceTaxId).filter((t) => t != null)
  );
  const allTaxIds = perContigOwnTaxIds.flat();
  if (allTaxIds.length === 0) {
    return { consensusTaxId: null, consensusRank: null, consensusName: null, perContigDistance: new Map() };
  }

  const consensusTaxId = tree.lca(allTaxIds);
  const consensusNode = consensusTaxId !== null ? tree.node(consensusTaxId) : null;

  const perContigDistance = new Map();
  binContigs.forEach((contig, i) => {
    const ownTaxIds = perContigOwnTaxIds[i];
    const otherTaxIds = perContigOwnTaxIds.filter((_, j) => j !== i).flat();
    if (ownTaxIds.length === 0 || otherTaxIds.length === 0) {
      perContigDistance.set(contig.id, null); // no signal without markers on both "this contig" and "the rest"
      return;
    }
    const othersLca = tree.lca(otherTaxIds);
    const othersNode = othersLca !== null ? tree.node(othersLca) : null;
    const combinedLca = tree.lca([...ownTaxIds, othersLca]);
    const combinedNode = combinedLca !== null ? tree.node(combinedLca) : null;
    perContigDistance.set(contig.id, othersNode && combinedNode ? othersNode.depth - combinedNode.depth : null);
  });

  return {
    consensusTaxId, consensusRank: consensusNode ? consensusNode.rank : null,
    consensusName: consensusNode ? consensusNode.name : null, perContigDistance,
  };
}

const exportsObj = { buildTaxonomyTreeFromLineage, loadTaxonomyTree, computeBinTaxonomicConsistency };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.markerTaxonomy = exportsObj;
}
})();
