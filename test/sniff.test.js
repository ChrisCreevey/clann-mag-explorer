const { test, report, assert } = require('./harness');
const { sniff } = require('../src/parsers/sniff');

test('sniffs a breport-shaped file', () => {
  const text = '100.00\t500\t500\tR\t1\troot\n50.00\t250\t250\tD\t2\tBacteria\n';
  assert.strictEqual(sniff(text).format, 'breport');
});

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
