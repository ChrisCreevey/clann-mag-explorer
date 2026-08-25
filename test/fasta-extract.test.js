const { test, report, assert } = require('./harness');
const { computeSequenceByteSpan, planFastaExtraction } = require('../src/model/fasta-extract');

test('computeSequenceByteSpan: single full line, no partial remainder', () => {
  // 60 bases per line, 1-byte line ending -> lineWidth 61; contig exactly 60 bases
  const span = computeSequenceByteSpan({ offset: 100, length: 60, lineBases: 60, lineWidth: 61, uniform: true });
  assert.deepStrictEqual(span, { offset: 100, byteLength: 61 });
});

test('computeSequenceByteSpan: multiple full lines plus a short last line', () => {
  // length 150 over lineBases 60 -> 2 full lines (120) + 30-base last line
  const span = computeSequenceByteSpan({ offset: 0, length: 150, lineBases: 60, lineWidth: 61, uniform: true });
  assert.deepStrictEqual(span, { offset: 0, byteLength: 153 }); // 61+61+31
});

test('computeSequenceByteSpan: exactly two full lines, no remainder', () => {
  const span = computeSequenceByteSpan({ offset: 0, length: 120, lineBases: 60, lineWidth: 61, uniform: true });
  assert.deepStrictEqual(span, { offset: 0, byteLength: 122 });
});

test('computeSequenceByteSpan: CRLF line endings (lineWidth - lineBases = 2)', () => {
  const span = computeSequenceByteSpan({ offset: 0, length: 70, lineBases: 60, lineWidth: 62, uniform: true });
  // 1 full line (60 bases, 62 bytes) + 10-base last line + 2-byte ending
  assert.deepStrictEqual(span, { offset: 0, byteLength: 74 });
});

test('computeSequenceByteSpan returns null for non-uniform wrapping', () => {
  const span = computeSequenceByteSpan({ offset: 0, length: 100, lineBases: 60, lineWidth: 61, uniform: false });
  assert.strictEqual(span, null);
});

test('computeSequenceByteSpan returns null for a zero-length contig', () => {
  const span = computeSequenceByteSpan({ offset: 0, length: 0, lineBases: 60, lineWidth: 61, uniform: true });
  assert.strictEqual(span, null);
});

test('planFastaExtraction groups entries per group, preserving contig order, and reports skips', () => {
  const records = [
    { id: 'c1', header: 'c1 desc', faiEntry: { offset: 0, length: 60, lineBases: 60, lineWidth: 61, uniform: true } },
    { id: 'c2', header: 'c2', faiEntry: { offset: 61, length: 100, lineBases: 60, lineWidth: 61, uniform: false } },
    { id: 'c3', header: 'c3', faiEntry: { offset: 200, length: 30, lineBases: 60, lineWidth: 61, uniform: true } },
  ];
  const contigIdsByGroup = new Map([
    ['MAG_1', ['c1', 'c2']],
    ['MAG_2', ['c3', 'ghost']],
  ]);
  const plan = planFastaExtraction(records, contigIdsByGroup);

  const mag1 = plan.get('MAG_1');
  assert.strictEqual(mag1.entries.length, 1);
  assert.strictEqual(mag1.entries[0].id, 'c1');
  assert.deepStrictEqual(mag1.skippedContigIds, ['c2']); // non-uniform, skipped

  const mag2 = plan.get('MAG_2');
  assert.strictEqual(mag2.entries.length, 1);
  assert.strictEqual(mag2.entries[0].id, 'c3');
  assert.deepStrictEqual(mag2.skippedContigIds, ['ghost']); // not in recordsById
});

report();
