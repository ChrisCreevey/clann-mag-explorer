const { test, report, assert } = require('./harness');
const { parseContigBinTable, looksLikeContigBinTable } = require('../src/parsers/contig-bin-table');

test('parses a DAS_Tool-style tab-delimited table with no header', () => {
  const text = 'contig_1\tbin.1\ncontig_2\tbin.1\ncontig_3\tbin.2\n';
  const rows = parseContigBinTable(text);
  assert.deepStrictEqual(rows, [
    { contigId: 'contig_1', binId: 'bin.1' },
    { contigId: 'contig_2', binId: 'bin.1' },
    { contigId: 'contig_3', binId: 'bin.2' },
  ]);
});

test('parses a CONCOCT-style comma-delimited table with a header row', () => {
  const text = 'contig_id,cluster_id\ncontig_1,0\ncontig_2,0\ncontig_3,1\n';
  const rows = parseContigBinTable(text);
  assert.deepStrictEqual(rows, [
    { contigId: 'contig_1', binId: '0' },
    { contigId: 'contig_2', binId: '0' },
    { contigId: 'contig_3', binId: '1' },
  ]);
});

test('skips blank lines and rows with fewer than 2 fields', () => {
  const text = 'contig_1\tbin.1\n\ncontig_2\n\ncontig_3\tbin.2\n';
  const rows = parseContigBinTable(text);
  assert.strictEqual(rows.length, 2);
  assert.strictEqual(rows[0].contigId, 'contig_1');
  assert.strictEqual(rows[1].contigId, 'contig_3');
});

test('handles Windows line endings', () => {
  const text = 'contig_1\tbin.1\r\ncontig_2\tbin.2\r\n';
  const rows = parseContigBinTable(text);
  assert.strictEqual(rows.length, 2);
  assert.strictEqual(rows[1].binId, 'bin.2');
});

test('empty input returns an empty array', () => {
  assert.deepStrictEqual(parseContigBinTable(''), []);
});

test('looksLikeContigBinTable accepts a valid table, headered or not', () => {
  assert.ok(looksLikeContigBinTable('contig_1\tbin.1\ncontig_2\tbin.2\n'));
  assert.ok(looksLikeContigBinTable('contig_id,cluster_id\ncontig_1,0\ncontig_2,1\n'));
});

test('looksLikeContigBinTable rejects a 6-column breport-shaped file', () => {
  const breportLine = '50.00\t100\t100\tS\t12345\ttest species\n';
  assert.strictEqual(looksLikeContigBinTable(breportLine.repeat(3)), false);
});

test('looksLikeContigBinTable rejects an inconsistent row shape', () => {
  assert.strictEqual(looksLikeContigBinTable('contig_1\tbin.1\ncontig_2\tbin.1\textra\n'), false);
});

report();
