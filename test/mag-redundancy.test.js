const { test, report, assert } = require('./harness');
const { computeMagRedundancy } = require('../src/model/mag-redundancy');

function contig(id, gcContent, composition) {
  return { id, gcContent, composition };
}

test('computeMagRedundancy flags two MAGs with near-identical composition as likely duplicate', () => {
  const records = [
    contig('a1', 0.5, { AAAA: 10, TTTT: 10 }),
    contig('a2', 0.5, { AAAA: 12, TTTT: 8 }),
    contig('b1', 0.5, { AAAA: 10, TTTT: 10 }),
    contig('b2', 0.5, { AAAA: 12, TTTT: 8 }),
  ];
  const putativeMags = [
    { magId: 'MAG_1', coreContigIds: ['a1', 'a2'], disputedContigIds: [] },
    { magId: 'MAG_2', coreContigIds: ['b1', 'b2'], disputedContigIds: [] },
  ];
  const pairs = computeMagRedundancy(records, putativeMags);
  assert.strictEqual(pairs.length, 1);
  assert.strictEqual(pairs[0].likelyDuplicate, true);
  assert.ok(pairs[0].compositionSimilarity > 0.99);
});

test('computeMagRedundancy does not flag MAGs with clearly different composition', () => {
  const records = [
    contig('a1', 0.3, { AAAA: 20, TTTT: 1 }),
    contig('a2', 0.3, { AAAA: 18, TTTT: 2 }),
    contig('b1', 0.7, { AAAA: 1, TTTT: 20 }),
    contig('b2', 0.7, { AAAA: 2, TTTT: 18 }),
  ];
  const putativeMags = [
    { magId: 'MAG_1', coreContigIds: ['a1', 'a2'], disputedContigIds: [] },
    { magId: 'MAG_2', coreContigIds: ['b1', 'b2'], disputedContigIds: [] },
  ];
  const pairs = computeMagRedundancy(records, putativeMags);
  assert.strictEqual(pairs.length, 1);
  assert.strictEqual(pairs[0].likelyDuplicate, false);
});

test('computeMagRedundancy skips MAGs below minContigs', () => {
  const records = [
    contig('a1', 0.5, { AAAA: 10 }),
    contig('b1', 0.5, { AAAA: 10 }),
    contig('b2', 0.5, { AAAA: 10 }),
  ];
  const putativeMags = [
    { magId: 'MAG_1', coreContigIds: ['a1'], disputedContigIds: [] }, // only 1 contig
    { magId: 'MAG_2', coreContigIds: ['b1', 'b2'], disputedContigIds: [] },
  ];
  const pairs = computeMagRedundancy(records, putativeMags);
  assert.strictEqual(pairs.length, 0);
});

test('computeMagRedundancy ranks pairs most-similar first', () => {
  const records = [
    contig('a1', 0.5, { AAAA: 10, TTTT: 10 }),
    contig('a2', 0.5, { AAAA: 10, TTTT: 10 }),
    contig('b1', 0.5, { AAAA: 10, TTTT: 10 }),
    contig('b2', 0.5, { AAAA: 10, TTTT: 10 }),
    contig('c1', 0.5, { AAAA: 1, TTTT: 30 }),
    contig('c2', 0.5, { AAAA: 1, TTTT: 30 }),
  ];
  const putativeMags = [
    { magId: 'MAG_1', coreContigIds: ['a1', 'a2'], disputedContigIds: [] },
    { magId: 'MAG_2', coreContigIds: ['b1', 'b2'], disputedContigIds: [] },
    { magId: 'MAG_3', coreContigIds: ['c1', 'c2'], disputedContigIds: [] },
  ];
  const pairs = computeMagRedundancy(records, putativeMags);
  assert.strictEqual(pairs.length, 3);
  assert.strictEqual(pairs[0].magIdA, 'MAG_1');
  assert.strictEqual(pairs[0].magIdB, 'MAG_2');
  assert.ok(pairs[0].compositionSimilarity >= pairs[1].compositionSimilarity);
  assert.ok(pairs[1].compositionSimilarity >= pairs[2].compositionSimilarity);
});

report();
