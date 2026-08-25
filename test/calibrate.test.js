const { test, report, assert } = require('./harness');
const { evenSample, percentile } = require('../build/03-calibrate');

test('evenSample returns the array unchanged when it is already <= n', () => {
  assert.deepStrictEqual(evenSample([1, 2, 3], 5), [1, 2, 3]);
});

test('evenSample picks exactly n evenly-spaced elements, spanning the full array', () => {
  const arr = Array.from({ length: 100 }, (_, i) => i);
  const sample = evenSample(arr, 10);
  assert.strictEqual(sample.length, 10);
  assert.strictEqual(sample[0], 0);
  // Deterministic — same input always yields the same sample (build reproducibility).
  assert.deepStrictEqual(sample, evenSample(arr, 10));
});

test('percentile returns null for an empty array', () => {
  assert.strictEqual(percentile([], 0.1), null);
});

test('percentile(sorted, 0) returns the minimum, matching a p10-style "safe floor" reading', () => {
  const sorted = [5, 10, 15, 20, 100];
  assert.strictEqual(percentile(sorted, 0), 5);
});

test('percentile respects sort order — caller must pre-sort ascending', () => {
  const sorted = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
  assert.strictEqual(percentile(sorted, 0.1), sorted[1]); // floor(0.1*10) = 1
});

report();
