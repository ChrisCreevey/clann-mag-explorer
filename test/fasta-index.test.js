const zlib = require('zlib');
const { test, report, assert } = require('./harness');
const { streamFasta, isGzipped } = require('../src/model/fasta-index');

const SAMPLE = [
  '>contig1 some description',
  'ACGTACGTAC',
  'GTACGTACGT',
  'ACGT',
  '>contig2',
  'GGGGCCCCGGGGCCCC',
  '',
].join('\n');

test('streamFasta parses multiple records in order with correct ids and lengths', async () => {
  const blob = new Blob([SAMPLE]);
  const records = [];
  const summary = await streamFasta(blob, (r) => records.push(r));

  assert.strictEqual(records.length, 2);
  assert.strictEqual(records[0].id, 'contig1');
  assert.strictEqual(records[0].header, 'contig1 some description');
  assert.strictEqual(records[0].length, 24);
  assert.strictEqual(records[1].id, 'contig2');
  assert.strictEqual(records[1].length, 16);
  assert.strictEqual(summary.contigCount, 2);
  assert.strictEqual(summary.totalLength, 40);
  assert.strictEqual(summary.sourceCompressed, false);
});

test('streamFasta builds a uniform .fai-style entry (lineBases/lineWidth/offset)', async () => {
  const blob = new Blob([SAMPLE]);
  const records = [];
  await streamFasta(blob, (r) => records.push(r));

  const c1 = records[0].faiEntry;
  // '>contig1 some description\n' is 26 bytes, so the sequence starts at byte 26.
  assert.strictEqual(c1.offset, 26);
  assert.strictEqual(c1.lineBases, 10);
  assert.strictEqual(c1.lineWidth, 11); // 10 bases + \n
  assert.strictEqual(c1.uniform, true); // last line (4 bases) shorter than the rest is still uniform
});

test('streamFasta flags a non-uniform line layout', async () => {
  const irregular = ['>c1', 'ACGT', 'AC', 'ACGTACGT', ''].join('\n');
  const blob = new Blob([irregular]);
  const records = [];
  await streamFasta(blob, (r) => records.push(r));
  assert.strictEqual(records[0].faiEntry.uniform, false);
});

test('streamFasta processes each contig with real per-contig stats attached', async () => {
  const blob = new Blob([SAMPLE]);
  const records = [];
  await streamFasta(blob, (r) => records.push(r));
  assert.strictEqual(records[1].gcContent, 1); // contig2 is all G/C
  assert.ok(records[0].composition && typeof records[0].composition === 'object');
  assert.ok(records[0].codingDensity >= 0 && records[0].codingDensity <= 1);
});

test('isGzipped detects gzip magic bytes', async () => {
  const gz = zlib.gzipSync(Buffer.from(SAMPLE));
  const plainBlob = new Blob([SAMPLE]);
  const gzBlob = new Blob([gz]);
  assert.strictEqual(await isGzipped(plainBlob), false);
  assert.strictEqual(await isGzipped(gzBlob), true);
});

test('streamFasta transparently decompresses a gzip-compressed FASTA', async () => {
  const gz = zlib.gzipSync(Buffer.from(SAMPLE));
  const gzBlob = new Blob([gz]);
  const records = [];
  const summary = await streamFasta(gzBlob, (r) => records.push(r));

  assert.strictEqual(records.length, 2);
  assert.strictEqual(records[0].id, 'contig1');
  assert.strictEqual(summary.sourceCompressed, true);
  assert.strictEqual(records[0].faiEntry.sourceCompressed, true);
});

test('streamFasta handles a file with no trailing newline', async () => {
  const noTrailingNewline = '>only\nACGTACGT';
  const blob = new Blob([noTrailingNewline]);
  const records = [];
  await streamFasta(blob, (r) => records.push(r));
  assert.strictEqual(records.length, 1);
  assert.strictEqual(records[0].length, 8);
});

report();
