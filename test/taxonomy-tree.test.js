const { test, report, assert } = require('./harness');
const { TaxonomyTree } = require('../src/model/taxonomy-tree');

test('getOrCreateNode builds parent-index chain', () => {
  const tree = new TaxonomyTree();
  const rootIdx = tree.getOrCreateNode(1, 'root', 'R', 0, null);
  const childIdx = tree.getOrCreateNode(2, 'child', 'D', 1, 1);
  assert.strictEqual(tree.parentIndex[childIdx], rootIdx);
  assert.strictEqual(tree.node(2).parentTaxid, 1);
});

test('lca finds common ancestor of a sibling pair', () => {
  const tree = new TaxonomyTree();
  tree.getOrCreateNode(1, 'root', 'R', 0, null);
  tree.getOrCreateNode(2, 'phylumA', 'P', 1, 1);
  tree.getOrCreateNode(3, 'genusA', 'G', 2, 2);
  tree.getOrCreateNode(4, 'genusB', 'G', 2, 2);
  assert.strictEqual(tree.lca([3, 4]), 2);
});

test('lca of identical taxid returns itself', () => {
  const tree = new TaxonomyTree();
  tree.getOrCreateNode(1, 'root', 'R', 0, null);
  tree.getOrCreateNode(2, 'genus', 'G', 1, 1);
  assert.strictEqual(tree.lca([2, 2]), 2);
});

test('lca across disjoint roots returns null', () => {
  const tree = new TaxonomyTree();
  tree.getOrCreateNode(1, 'rootA', 'R', 0, null);
  tree.getOrCreateNode(2, 'rootB', 'R', 0, null);
  assert.strictEqual(tree.lca([1, 2]), null);
});

test('lca of a single taxid returns itself', () => {
  const tree = new TaxonomyTree();
  tree.getOrCreateNode(1, 'root', 'R', 0, null);
  tree.getOrCreateNode(2, 'genus', 'G', 1, 1);
  assert.strictEqual(tree.lca([2]), 2);
});

report();
