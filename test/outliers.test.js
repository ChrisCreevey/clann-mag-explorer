const { test, report, assert } = require('./harness');
const { centroid, euclideanDistance, computeCentroidZScores, computeBinOutliers } = require('../src/model/outliers');

test('centroid is the mean of each dimension', () => {
  assert.deepStrictEqual(centroid([[0, 0], [2, 4]]), [1, 2]);
});

test('euclideanDistance of identical points is 0', () => {
  assert.strictEqual(euclideanDistance([1, 2, 3], [1, 2, 3]), 0);
});

test('euclideanDistance matches the 3-4-5 triangle', () => {
  assert.strictEqual(euclideanDistance([0, 0], [3, 4]), 5);
});

test('computeCentroidZScores gives an outlier a clearly higher z than tightly clustered points', () => {
  const vectorById = new Map([
    ['a', [0, 0]], ['b', [0.1, 0.1]], ['c', [-0.1, 0.1]], ['d', [10, 10]], // d is far from the others
  ]);
  const z = computeCentroidZScores(['a', 'b', 'c', 'd'], vectorById);
  assert.ok(z.get('d') > z.get('a'));
  assert.ok(z.get('d') > z.get('b'));
});

test('computeCentroidZScores returns 0 for every contig when all points are identical (zero spread)', () => {
  const vectorById = new Map([['a', [1, 1]], ['b', [1, 1]], ['c', [1, 1]]]);
  const z = computeCentroidZScores(['a', 'b', 'c'], vectorById);
  assert.strictEqual(z.get('a'), 0);
  assert.strictEqual(z.get('c'), 0);
});

test('computeBinOutliers flags a compositional outlier using only composition when no coverage is loaded', () => {
  const binContigs = [
    { id: 'a', composition: { AAAA: 0.5, CCCC: 0.5 } },
    { id: 'b', composition: { AAAA: 0.52, CCCC: 0.48 } },
    { id: 'c', composition: { AAAA: 0.01, CCCC: 0.99 } }, // the odd one out
  ];
  const result = computeBinOutliers(binContigs);
  assert.strictEqual(result.get('a').coverageZ, null);
  assert.ok(result.get('c').compositionZ > result.get('a').compositionZ);
  assert.strictEqual(result.get('c').combinedZ, result.get('c').compositionZ);
});

test('computeBinOutliers uses coverage too when every contig in the bin has it, combined as max(compositionZ, coverageZ)', () => {
  const binContigs = [
    { id: 'a', composition: { AAAA: 0.5, CCCC: 0.5 }, coverageDepths: [10, 10] },
    { id: 'b', composition: { AAAA: 0.5, CCCC: 0.5 }, coverageDepths: [10.1, 9.9] },
    { id: 'c', composition: { AAAA: 0.5, CCCC: 0.5 }, coverageDepths: [80, 80] }, // coverage outlier, composition-identical
  ];
  const result = computeBinOutliers(binContigs);
  assert.ok(result.get('c').coverageZ > result.get('a').coverageZ);
  assert.strictEqual(result.get('c').compositionZ, 0); // identical composition -> no compositional signal
  assert.strictEqual(result.get('c').combinedZ, result.get('c').coverageZ); // max() picked coverage here
});

test('computeBinOutliers falls back to composition-only if even one contig lacks coverage data', () => {
  const binContigs = [
    { id: 'a', composition: { AAAA: 0.5, CCCC: 0.5 }, coverageDepths: [10, 10] },
    { id: 'b', composition: { AAAA: 0.5, CCCC: 0.5 } }, // no coverageDepths at all
  ];
  const result = computeBinOutliers(binContigs);
  assert.strictEqual(result.get('a').coverageZ, null);
  assert.strictEqual(result.get('b').coverageZ, null);
});

report();
