const { test, report, assert } = require('./harness');
const { findOrfRanges, highlightOrfsHtml } = require('../src/model/orf-finder');
const { reverseComplement, translateFrame } = require('../src/model/translate');

test('findOrfRanges finds a forward-strand ORF between two stops, ignores it below minAaLength', () => {
  // TAA (stop) + 150x GGC (Gly) + TAA (stop): a clean 150-aa run in frame
  // 0. GGC repeated happens to contain no stop codon in frames 1/2
  // either (a property of this particular codon, not something worth
  // engineering around), so those show up as their own long ORFs too —
  // real, correct behaviour, not a bug — hence checking the frame-0
  // range *exists* rather than asserting it's the only result.
  const seq = 'TAA' + 'GGC'.repeat(150) + 'TAA';
  const ranges = findOrfRanges(seq, 100);
  const found = ranges.find((r) => r.strand === 'fwd' && r.start === 3 && r.end === 453);
  assert.ok(found, 'expected a forward ORF spanning the 150 Gly codons between the two stops');
  assert.strictEqual((found.end - found.start) / 3, 150);

  const tooStrict = findOrfRanges(seq, 200);
  assert.ok(!tooStrict.some((r) => r.strand === 'fwd'), 'a 150-aa run (and the ~151-aa runs in the other two frames) should not pass a 200-aa minimum');
});

test('findOrfRanges finds a reverse-strand ORF at the correct forward-sequence coordinates', () => {
  // Build the *forward* translation of the intended reverse-complement
  // strand directly (stop, 80x Gly, stop), then take its reverse
  // complement as the actual (forward-oriented) sequence under test —
  // the mirror-image of the forward-ORF test above.
  const rcStrand = 'TAA' + 'GGC'.repeat(80) + 'TAA';
  const seq = reverseComplement(rcStrand);

  const ranges = findOrfRanges(seq, 50);
  // Frame 0 (the one this test targets) produces an exact 80-aa/240nt
  // segment, bounded by the two real stop codons; frames 1/2 may also
  // independently qualify (no internal stop, same property GGC repeated
  // has in the forward test above) but won't share this exact length, so
  // searching for it is robust regardless of how many other reverse
  // ranges are also found.
  const found = ranges.find((r) => r.strand === 'rev' && r.end - r.start === 240);
  assert.ok(found, 'expected a reverse-strand ORF spanning exactly the 80 Gly codons between the two stops');

  // Verify it independently of how it was derived: slicing it out and
  // re-translating its reverse complement should be a clean 80-residue
  // run of Gly with no stops, regardless of the exact offset arithmetic
  // used internally.
  const segment = seq.slice(found.start, found.end);
  const aa = translateFrame(reverseComplement(segment), 0);
  assert.strictEqual(aa, 'G'.repeat(80));
});

test('findOrfRanges treats an ambiguous-base codon (X) as a segment boundary, same as a stop', () => {
  const seq = 'GGC'.repeat(100) + 'NNN' + 'GGC'.repeat(100); // 300nt Gly, one N-codon, 300nt Gly
  const ranges = findOrfRanges(seq, 100).filter((r) => r.strand === 'fwd');
  assert.strictEqual(ranges.length, 2, 'the ambiguous codon should split one long run into two ORFs, not one');
  assert.strictEqual((ranges[0].end - ranges[0].start) / 3, 100);
  assert.strictEqual((ranges[1].end - ranges[1].start) / 3, 100);
});

test('findOrfRanges returns nothing for a sequence with no run long enough on any frame', () => {
  const seq = 'TAAGGCTAAGGCTAAGGCTAA'; // stops every few codons, nothing close to 100aa
  assert.deepStrictEqual(findOrfRanges(seq, 100), []);
});

test('highlightOrfsHtml wraps a highlighted range in a span and still line-wraps at the given width', () => {
  const seq = 'A'.repeat(10) + 'C'.repeat(10) + 'A'.repeat(10); // 30nt, middle 10 highlighted
  const html = highlightOrfsHtml(seq, [{ start: 10, end: 20, strand: 'fwd' }], 12);

  // Round-trips back to the original sequence once markup/newlines are stripped.
  assert.strictEqual(html.replace(/<[^>]+>/g, '').replace(/\n/g, ''), seq);

  // Newlines land at fixed positions (every 12 visible characters) regardless of span boundaries.
  const withoutTags = html.replace(/<\/?span[^>]*>/g, '');
  assert.strictEqual(withoutTags, seq.slice(0, 12) + '\n' + seq.slice(12, 24) + '\n' + seq.slice(24));

  // The highlighted run is wrapped in exactly one orf-fwd span containing only the highlighted characters.
  const spanMatch = html.match(/<span class="orf-fwd">([\s\S]*?)<\/span>/);
  assert.ok(spanMatch, 'expected one orf-fwd span');
  assert.strictEqual(spanMatch[1].replace(/\n/g, ''), 'C'.repeat(10));
});

test('highlightOrfsHtml marks an overlapping forward+reverse region as orf-both', () => {
  const seq = 'A'.repeat(20);
  const html = highlightOrfsHtml(seq, [
    { start: 5, end: 15, strand: 'fwd' },
    { start: 10, end: 20, strand: 'rev' },
  ], 100);
  assert.ok(html.includes('<span class="orf-fwd">AAAAA</span>'), 'positions 5-9: forward only');
  assert.ok(html.includes('<span class="orf-both">AAAAA</span>'), 'positions 10-14: both strands');
  assert.ok(html.includes('<span class="orf-rev">AAAAA</span>'), 'positions 15-19: reverse only');
});

test('highlightOrfsHtml returns the plain sequence, still wrapped, when no ranges are given', () => {
  const seq = 'ACGT'.repeat(10); // 40nt
  assert.strictEqual(highlightOrfsHtml(seq, [], 20), 'ACGT'.repeat(5) + '\n' + 'ACGT'.repeat(5));
});

report();
