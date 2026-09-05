const { test, report, assert } = require('./harness');
const { sniff } = require('../src/parsers/sniff');

test('sniffs a DAS_Tool-style contig-bin table', () => {
  const text = 'contig_1\tbin.1\ncontig_2\tbin.2\n';
  assert.strictEqual(sniff(text).format, 'contig-bin-table');
});

test('sniffs a CONCOCT-style headered contig-bin table', () => {
  const text = 'contig_id,cluster_id\ncontig_1,0\ncontig_2,1\n';
  assert.strictEqual(sniff(text).format, 'contig-bin-table');
});

test('returns unknown for content matching neither shape', () => {
  const text = 'this is not tabular data at all\njust some prose\n';
  const result = sniff(text);
  assert.strictEqual(result.format, 'unknown');
  assert.ok(result.reason);
});

report();
