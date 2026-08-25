const { test, report, assert } = require('./harness');
const { buildTaxonomyTreeFromLineage, computeBinTaxonomicConsistency } = require('../src/model/marker-taxonomy');

// A small synthetic lineage: root(1) -> Bacteria(2) -> { Proteobacteria(3) -> E.coli(10), Firmicutes(4) -> Bacillus(11) }
const LINEAGE = [
  { taxid: 1, parentTaxid: 1, rank: 'R', name: 'root' },
  { taxid: 2, parentTaxid: 1, rank: 'D', name: 'Bacteria' },
  { taxid: 3, parentTaxid: 2, rank: 'P', name: 'Proteobacteria' },
  { taxid: 4, parentTaxid: 2, rank: 'P', name: 'Firmicutes' },
  { taxid: 10, parentTaxid: 3, rank: 'S', name: 'Escherichia coli' },
  { taxid: 11, parentTaxid: 4, rank: 'S', name: 'Bacillus subtilis' },
];

test('buildTaxonomyTreeFromLineage builds correct parent/depth regardless of input row order', () => {
  const shuffled = [LINEAGE[4], LINEAGE[1], LINEAGE[0], LINEAGE[5], LINEAGE[3], LINEAGE[2]];
  const tree = buildTaxonomyTreeFromLineage(shuffled);
  assert.strictEqual(tree.node(10).parentTaxid, 3);
  assert.strictEqual(tree.node(10).depth, 3);
  assert.strictEqual(tree.node(1).depth, 0);
});

test('computeBinTaxonomicConsistency: all markers from the same lineage -> consensus resolves at species, all contigs distance 0', () => {
  const tree = buildTaxonomyTreeFromLineage(LINEAGE);
  const binContigs = [
    { id: 'c1', markerHits: [{ provenanceTaxId: 10 }] },
    { id: 'c2', markerHits: [{ provenanceTaxId: 10 }] },
  ];
  const result = computeBinTaxonomicConsistency(binContigs, tree);
  assert.strictEqual(result.consensusTaxId, 10);
  assert.strictEqual(result.consensusRank, 'S');
  assert.strictEqual(result.perContigDistance.get('c1'), 0);
  assert.strictEqual(result.perContigDistance.get('c2'), 0);
});

test('computeBinTaxonomicConsistency: a chimeric bin flags the odd contig out via leave-one-out distance', () => {
  const tree = buildTaxonomyTreeFromLineage(LINEAGE);
  const binContigs = [
    { id: 'c1', markerHits: [{ provenanceTaxId: 10 }] }, // E. coli
    { id: 'c2', markerHits: [{ provenanceTaxId: 10 }] }, // E. coli
    { id: 'c3', markerHits: [{ provenanceTaxId: 11 }] }, // Bacillus — the odd one out
  ];
  const result = computeBinTaxonomicConsistency(binContigs, tree);
  // whole-bin consensus = LCA(10,10,11) = Bacteria (taxid 2)
  assert.strictEqual(result.consensusTaxId, 2);
  assert.strictEqual(result.consensusRank, 'D');
  // c1: the other two contigs are {c2=10, c3=11} -> their own LCA is Bacteria(2, depth 1);
  // adding c1's own 10 back in doesn't pull it any further (10 is already under 2) -> distance 0
  assert.strictEqual(result.perContigDistance.get('c1'), 0);
  // c3: the other two contigs are {c1=10, c2=10} -> their own LCA is E. coli (10, depth 3);
  // adding c3's own 11 pulls that up to Bacteria (2, depth 1) -> distance 3-1 = 2, clearly flagged
  assert.strictEqual(result.perContigDistance.get('c3'), 2);
});

test('computeBinTaxonomicConsistency: a contig with no marker hits gets a null distance, not zero', () => {
  const tree = buildTaxonomyTreeFromLineage(LINEAGE);
  const binContigs = [
    { id: 'c1', markerHits: [{ provenanceTaxId: 10 }] },
    { id: 'c2', markerHits: [] },
  ];
  const result = computeBinTaxonomicConsistency(binContigs, tree);
  assert.strictEqual(result.perContigDistance.get('c2'), null);
});

test('computeBinTaxonomicConsistency: a bin with no marker hits at all returns a null consensus', () => {
  const tree = buildTaxonomyTreeFromLineage(LINEAGE);
  const result = computeBinTaxonomicConsistency([{ id: 'c1', markerHits: [] }], tree);
  assert.strictEqual(result.consensusTaxId, null);
});

report();
