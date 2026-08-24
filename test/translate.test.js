const { test, report, assert } = require('./harness');
const { reverseComplement, translateFrame, translateSixFrames } = require('../src/model/translate');

test('reverseComplement reverses and complements', () => {
  assert.strictEqual(reverseComplement('ATGC'), 'GCAT');
  assert.strictEqual(reverseComplement('AAAA'), 'TTTT');
});

test('reverseComplement maps ambiguity codes to N', () => {
  assert.strictEqual(reverseComplement('ATGN'), 'NCAT');
});

test('translateFrame reads codons from the given offset', () => {
  // ATG GCA TTT -> M A F
  assert.strictEqual(translateFrame('ATGGCATTT', 0), 'MAF');
  // shifted by 1: TGG CAT TT(partial dropped)
  assert.strictEqual(translateFrame('ATGGCATTT', 1), 'WH');
});

test('translateFrame marks stop codons with a sentinel', () => {
  // ATG TAA GGG -> M * G
  assert.strictEqual(translateFrame('ATGTAAGGG', 0), 'M*G');
});

test('translateFrame marks non-ACGT codons as X, not a stop', () => {
  assert.strictEqual(translateFrame('ATGNNNGGG', 0), 'MXG');
});

test('translateSixFrames returns 3 forward + 3 reverse-complement frames', () => {
  const frames = translateSixFrames('ATGGCATTTATGGCATTTA');
  assert.strictEqual(frames.length, 6);
  assert.strictEqual(frames[0], translateFrame('ATGGCATTTATGGCATTTA', 0));
  assert.strictEqual(frames[3], translateFrame(reverseComplement('ATGGCATTTATGGCATTTA'), 0));
});

report();
