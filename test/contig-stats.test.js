const { test, report, assert } = require('./harness');
const {
  computeContigStats,
  computeTetranucleotideComposition,
  computeCodingDensity,
  getCanonicalKmerIndex,
} = require('../src/model/contig-stats');

test('computeContigStats reports length, GC content, GC skew', () => {
  const stats = computeContigStats('GGGGCCCC'); // 4 G, 4 C
  assert.strictEqual(stats.length, 8);
  assert.strictEqual(stats.gcContent, 1);
  assert.strictEqual(stats.gcSkew, 0); // (G-C)/(G+C) = 0
});

test('gcSkew is positive when G outnumbers C', () => {
  const stats = computeContigStats('GGGGCCCCAT'.replace('CCCC', 'CC')); // GGGGCCAT: 4G 2C
  assert.ok(stats.gcSkew > 0);
});

test('handles an all-N sequence without dividing by zero', () => {
  const stats = computeContigStats('NNNN');
  assert.strictEqual(stats.gcContent, 0);
  assert.strictEqual(stats.gcSkew, 0);
  assert.strictEqual(stats.ambiguousBaseCount, 4);
});

test('canonical k-mer index collapses 256 4-mers to 136 canonical bins', () => {
  const { canonicalList } = getCanonicalKmerIndex();
  assert.strictEqual(canonicalList.length, 136);
});

test('composition collapses a k-mer and its reverse complement to the same bin', () => {
  // AAAA and its revcomp TTTT must land in the same canonical bin.
  const compAAAA = computeTetranucleotideComposition('AAAA');
  const compTTTT = computeTetranucleotideComposition('TTTT');
  assert.deepStrictEqual(compAAAA, compTTTT);
});

test('composition frequencies sum to 1 when at least one full window exists', () => {
  const comp = computeTetranucleotideComposition('ACGTACGTAC');
  const sum = Object.values(comp).reduce((a, b) => a + b, 0);
  assert.ok(Math.abs(sum - 1) < 1e-9);
});

test('composition skips windows touching a non-ACGT base', () => {
  const comp = computeTetranucleotideComposition('ACGN');
  const sum = Object.values(comp).reduce((a, b) => a + b, 0);
  assert.strictEqual(sum, 0); // the only possible window (ACGN) is skipped
});

test('computeCodingDensity counts only segments at/above the minimum length', () => {
  // frame with one long run of 25 non-stop residues and one short run of 5
  const frame = 'X'.repeat(25) + '*' + 'X'.repeat(5);
  const density = computeCodingDensity([frame], 20);
  assert.strictEqual(density, 25 / frame.replace(/\*/g, '').length);
});

report();
