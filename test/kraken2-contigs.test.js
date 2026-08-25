const { test, report, assert } = require('./harness');
const { parseKraken2ContigCalls, looksLikeKraken2ContigOutput } = require('../src/parsers/kraken2-contigs');

test('parses bare-taxid classified and unclassified rows', () => {
  const text = 'C\tcontig_1\t562\t5000\t562:100\nU\tcontig_2\t0\t3000\t0:80\n';
  const calls = parseKraken2ContigCalls(text);
  assert.deepStrictEqual(calls, [
    { contigId: 'contig_1', classified: true, taxId: 562 },
    { contigId: 'contig_2', classified: false, taxId: 0 },
  ]);
});

test('parses --use-names format ("name (taxid N)")', () => {
  const text = 'C\tcontig_1\tEscherichia coli (taxid 562)\t5000\t562:100\n';
  const calls = parseKraken2ContigCalls(text);
  assert.strictEqual(calls[0].taxId, 562);
});

test('skips malformed lines rather than throwing', () => {
  const text = 'C\tcontig_1\t562\t5000\t562:100\nnot a valid line\nC\tcontig_2\t9\t100\t9:5\n';
  const calls = parseKraken2ContigCalls(text);
  assert.strictEqual(calls.length, 2);
});

test('looksLikeKraken2ContigOutput accepts a real-shaped file', () => {
  assert.ok(looksLikeKraken2ContigOutput('C\tcontig_1\t562\t5000\t562:100\nU\tcontig_2\t0\t3000\t0:80\n'));
});

test('looksLikeKraken2ContigOutput rejects a breport-shaped file', () => {
  const breportLine = '50.00\t100\t100\tS\t12345\ttest species\n';
  assert.strictEqual(looksLikeKraken2ContigOutput(breportLine.repeat(3)), false);
});

test('looksLikeKraken2ContigOutput rejects a contig-bin table', () => {
  assert.strictEqual(looksLikeKraken2ContigOutput('contig_1\tbin.1\ncontig_2\tbin.2\n'), false);
});

report();
