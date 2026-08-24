const { test, report, assert } = require('./harness');
const { reverseComplement, translateFrame, translateReverseFrame, translateSixFrames } = require('../src/model/translate');

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

test('translateReverseFrame matches translateFrame(reverseComplement(seq), offset) exactly', () => {
  // Deterministic pseudo-random sequences at several lengths, including
  // ones not evenly divisible by 3 and lengths shorter than one codon,
  // to exercise the boundary/off-by-one cases in the direct-lookup path.
  let seed = 42;
  function rand() { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; }
  const bases = 'ACGT';
  for (const length of [0, 1, 2, 3, 4, 5, 6, 7, 50, 301, 1000]) {
    let seq = '';
    for (let i = 0; i < length; i++) seq += bases[Math.floor(rand() * 4)];
    for (const offset of [0, 1, 2]) {
      const viaRC = translateFrame(reverseComplement(seq), offset);
      const direct = translateReverseFrame(seq, offset);
      assert.strictEqual(direct, viaRC, `mismatch at length=${length} offset=${offset}`);
    }
  }
});

test('translateReverseFrame treats ambiguous bases as X, same as the RC-string path', () => {
  const seq = 'ATGNCGTAA'; // includes an ambiguity code mid-sequence
  for (const offset of [0, 1, 2]) {
    assert.strictEqual(translateReverseFrame(seq, offset), translateFrame(reverseComplement(seq), offset));
  }
});

test('translateSixFrames reverse frames still match the reverse-complement definition', () => {
  const seq = 'ATGGCATTTATGGCATTTA';
  const frames = translateSixFrames(seq);
  const rc = reverseComplement(seq);
  assert.strictEqual(frames[3], translateFrame(rc, 0));
  assert.strictEqual(frames[4], translateFrame(rc, 1));
  assert.strictEqual(frames[5], translateFrame(rc, 2));
});

report();
