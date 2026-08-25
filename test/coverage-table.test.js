const { test, report, assert } = require('./harness');
const { parseCoverageTable, looksLikeCoverageTable } = require('../src/parsers/coverage-table');

test('parses a generic headerless table (contig + N depth columns)', () => {
  const text = 'contig_1\t12.5\t3.1\ncontig_2\t0.4\t9.9\n';
  const { sampleNames, rows } = parseCoverageTable(text);
  assert.deepStrictEqual(sampleNames, ['sample1', 'sample2']);
  assert.strictEqual(rows.length, 2);
  assert.deepStrictEqual(rows[0], { contigId: 'contig_1', depths: [12.5, 3.1] });
});

test('parses a generic headered table, using header cells as sample names', () => {
  const text = 'contig\tsiteA\tsiteB\ncontig_1\t12.5\t3.1\n';
  const { sampleNames, rows } = parseCoverageTable(text);
  assert.deepStrictEqual(sampleNames, ['siteA', 'siteB']);
  assert.strictEqual(rows[0].depths[0], 12.5);
});

test('parses MetaBAT2 jgi_summarize_bam_contig_depths output, dropping -var and derived columns', () => {
  const text = 'contigName\tcontigLen\ttotalAvgDepth\tsample1.bam\tsample1.bam-var\tsample2.bam\tsample2.bam-var\n' +
    'contig_1\t5000\t10.2\t12.5\t1.1\t8.0\t0.9\n';
  const { sampleNames, rows } = parseCoverageTable(text);
  assert.deepStrictEqual(sampleNames, ['sample1.bam', 'sample2.bam']);
  assert.deepStrictEqual(rows[0].depths, [12.5, 8.0]);
});

test('empty input returns no rows', () => {
  assert.deepStrictEqual(parseCoverageTable(''), { sampleNames: [], rows: [] });
});

test('looksLikeCoverageTable accepts a numeric multi-column table', () => {
  assert.ok(looksLikeCoverageTable('contig_1\t12.5\t3.1\ncontig_2\t0.4\t9.9\n'));
});

test('looksLikeCoverageTable rejects a table with non-numeric data columns', () => {
  assert.strictEqual(looksLikeCoverageTable('contig_1\tbin.1\ncontig_2\tbin.2\n'), false);
});

report();
