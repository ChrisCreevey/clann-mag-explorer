const { test, report, assert } = require('./harness');
const {
  parseIndexBinary, parseRefSeqsBinary, buildKeyLookup,
  extendUngapped, searchContigForMarkers,
} = require('../src/model/marker-genes');
const { forEachReducedKmer } = require('../src/model/reduced-alphabet');
const fixtures = require('./fixtures/marker-gene-sequences.json');

const K = 5;

// Builds an in-memory assets object with the same shape parseAssets()
// produces, straight from {family: [{header,taxId,seq}]} — the small-scale
// equivalent of build/01-cluster.js + build/02-index.js, uncapped (test
// data is tiny, no need for the real pipeline's per-key hit cap).
function buildTestAssets(byFamily) {
  const familyNames = Object.keys(byFamily).sort();
  const records = [];
  for (const family of familyNames) {
    for (const r of byFamily[family]) records.push({ family, ...r });
  }

  const seqOffsets = new Uint32Array(records.length + 1);
  const taxId = new Uint32Array(records.length);
  const familyIndex = new Uint8Array(records.length);
  let totalLen = 0;
  for (const r of records) totalLen += r.seq.length;
  const residues = new Uint8Array(totalLen);

  let cursor = 0;
  const familyIdxByName = new Map(familyNames.map((f, i) => [f, i]));
  for (let i = 0; i < records.length; i++) {
    const r = records[i];
    seqOffsets[i] = cursor;
    for (let j = 0; j < r.seq.length; j++) residues[cursor + j] = r.seq.charCodeAt(j);
    cursor += r.seq.length;
    taxId[i] = r.taxId;
    familyIndex[i] = familyIdxByName.get(r.family);
  }
  seqOffsets[records.length] = cursor;

  const hitsByKey = new Map();
  for (let seqId = 0; seqId < records.length; seqId++) {
    forEachReducedKmer(records[seqId].seq, K, (code, pos) => {
      if (!hitsByKey.has(code)) hitsByKey.set(code, []);
      hitsByKey.get(code).push([seqId, pos]);
    });
  }
  // Sparse (populated-keys-only) shape — matches parseIndexBinary's output
  // (see src/model/marker-genes.js), not a dense ALPHABET_SIZE^K array.
  const sortedCodes = [...hitsByKey.keys()].sort((a, b) => a - b);
  const sortedKeys = Uint32Array.from(sortedCodes);
  const keyOffsets = new Uint32Array(sortedCodes.length + 1);
  let hitCursor = 0;
  const hitRefSeqId = [];
  const hitPosition = [];
  for (let i = 0; i < sortedCodes.length; i++) {
    keyOffsets[i] = hitCursor;
    for (const [seqId, pos] of hitsByKey.get(sortedCodes[i])) { hitRefSeqId.push(seqId); hitPosition.push(pos); hitCursor++; }
  }
  keyOffsets[sortedCodes.length] = hitCursor;
  const keyLookup = buildKeyLookup(sortedKeys);

  return {
    index: {
      k: K, numPopulatedKeys: sortedCodes.length, numHits: hitCursor,
      sortedKeys, keyLookup, keyOffsets,
      hitRefSeqId: Uint16Array.from(hitRefSeqId), hitPosition: Uint16Array.from(hitPosition),
    },
    refSeqs: { numSeqs: records.length, seqOffsets, taxId, familyIndex, residues },
    familyNames,
  };
}

function padFrames(realFrame, others = 5) {
  // Six frames total; only frame 0 carries the real content, the rest are
  // short unrelated padding so a real search has something to iterate.
  const frames = [realFrame];
  for (let i = 0; i < others; i++) frames.push('MKLVAGTS');
  return frames;
}

// ---- Binary round-trip ----

test('parseIndexBinary/parseRefSeqsBinary round-trip a hand-built buffer', () => {
  // A tiny, deliberately simple sparse index: 1 populated key (code 2), 2 hits; 1 ref seq.
  const numPopulatedKeys = 1, numHits = 2, numSeqs = 1, residueBytes = 4;
  const idxHeader = Buffer.alloc(16);
  idxHeader.write('SCGI', 0, 'ascii');
  idxHeader.writeUInt8(2, 4); idxHeader.writeUInt8(K, 5); idxHeader.writeUInt8(10, 6); idxHeader.writeUInt8(2, 7);
  idxHeader.writeUInt32LE(numPopulatedKeys, 8); idxHeader.writeUInt32LE(numHits, 12);
  const sortedKeys = Uint32Array.from([2]);
  const keyOffsets = Uint32Array.from([0, 2]);
  const hitRefSeqId = Uint16Array.from([0, 0]);
  const hitPosition = Uint16Array.from([1, 3]);
  const idxBuf = Buffer.concat([idxHeader, Buffer.from(sortedKeys.buffer), Buffer.from(keyOffsets.buffer), Buffer.from(hitRefSeqId.buffer), Buffer.from(hitPosition.buffer)]);

  const parsed = parseIndexBinary(idxBuf.buffer.slice(idxBuf.byteOffset, idxBuf.byteOffset + idxBuf.length));
  assert.strictEqual(parsed.k, K);
  assert.strictEqual(parsed.numHits, numHits);
  assert.deepStrictEqual([...parsed.sortedKeys], [...sortedKeys]);
  assert.deepStrictEqual([...parsed.keyOffsets], [...keyOffsets]);
  assert.deepStrictEqual([...parsed.hitRefSeqId], [0, 0]);
  assert.deepStrictEqual([...parsed.hitPosition], [1, 3]);

  const refHeader = Buffer.alloc(16);
  refHeader.write('SCGR', 0, 'ascii');
  refHeader.writeUInt8(1, 4);
  refHeader.writeUInt32LE(numSeqs, 8); refHeader.writeUInt32LE(residueBytes, 12);
  const seqOffsets = Uint32Array.from([0, 4]);
  const taxId = Uint32Array.from([42]);
  const familyIndex = Uint8Array.from([0]);
  const residues = Uint8Array.from(Buffer.from('MSLK', 'ascii'));
  const refBuf = Buffer.concat([refHeader, Buffer.from(seqOffsets.buffer), Buffer.from(taxId.buffer), Buffer.from(familyIndex.buffer), Buffer.from(residues.buffer)]);

  const parsedRef = parseRefSeqsBinary(refBuf.buffer.slice(refBuf.byteOffset, refBuf.byteOffset + refBuf.length));
  assert.strictEqual(parsedRef.numSeqs, 1);
  assert.deepStrictEqual([...parsedRef.taxId], [42]);
  assert.strictEqual(Buffer.from(parsedRef.residues).toString('ascii'), 'MSLK');
});

test('parseIndexBinary rejects a bad magic number', () => {
  const buf = Buffer.alloc(16);
  buf.write('NOPE', 0, 'ascii');
  assert.throws(() => parseIndexBinary(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
});

// ---- extendUngapped ----

test('extendUngapped scores a perfect match at its full length', () => {
  const seq = 'MSLKCGIVGL';
  const residues = Uint8Array.from(Buffer.from(seq, 'ascii'));
  const result = extendUngapped(seq, residues, 0, seq.length, 5, 5, 15);
  assert.strictEqual(result.alignedLength, seq.length);
  assert.ok(result.score > 0);
});

test('extendUngapped halts at a stop-codon sentinel without special-case logic', () => {
  const query = 'AAAAA*BBBBB'; // '*' should stop extension from either side
  const ref =   'AAAAAXBBBBB';
  const residues = Uint8Array.from(Buffer.from(ref, 'ascii'));
  const result = extendUngapped(query, residues, 0, ref.length, 2, 2, 15);
  // Anchored inside the left block, extension should not cross the '*' at index 5
  assert.ok(result.frameEnd < 5, `expected extension to stop before the sentinel, got frameEnd=${result.frameEnd}`);
});

// ---- searchContigForMarkers: realistic fixtures ----

test('a contig containing a real COG0012 sequence is called as COG0012, not COG0016', () => {
  const assets = buildTestAssets(fixtures);
  const embedded = fixtures.cog0012[0].seq; // exact copy of one representative
  const frames = padFrames(embedded);

  const calls = searchContigForMarkers(frames, assets, { minRepresentatives: 2 });
  const families = calls.map((c) => c.family);
  assert.ok(families.includes('cog0012'), `expected cog0012 in ${JSON.stringify(families)}`);
  assert.ok(!families.includes('cog0016'), `did not expect cog0016 in ${JSON.stringify(families)}`);

  const call = calls.find((c) => c.family === 'cog0012');
  assert.ok(call.representativeCount >= 2, `expected multi-representative agreement, got ${call.representativeCount}`);
  assert.strictEqual(call.provenanceTaxId, fixtures.cog0012[0].taxId);
});

test('an unrelated random sequence is not called as any family', () => {
  const assets = buildTestAssets(fixtures);
  let seed = 7;
  function rand() { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; }
  const AA = 'ACDEFGHIKLMNPQRSTVWY';
  let random = '';
  for (let i = 0; i < 300; i++) random += AA[Math.floor(rand() * AA.length)];

  const calls = searchContigForMarkers(padFrames(random), assets);
  assert.strictEqual(calls.length, 0, `expected no calls, got ${JSON.stringify(calls)}`);
});

test('multi-representative agreement check rejects a family with only one qualifying representative', () => {
  const assets = buildTestAssets({ cog0012: [fixtures.cog0012[0]] }); // only 1 representative in the index
  const frames = padFrames(fixtures.cog0012[0].seq);
  const calls = searchContigForMarkers(frames, assets, { minRepresentatives: 2 });
  assert.strictEqual(calls.length, 0, 'a single representative should fail the multi-representative agreement check');
});

test('lowering minRepresentatives to 1 allows a single-representative call through', () => {
  const assets = buildTestAssets({ cog0012: [fixtures.cog0012[0]] });
  const frames = padFrames(fixtures.cog0012[0].seq);
  const calls = searchContigForMarkers(frames, assets, { minRepresentatives: 1 });
  assert.strictEqual(calls.length, 1);
  assert.strictEqual(calls[0].family, 'cog0012');
});

// ---- Per-family threshold overrides (build/03-calibrate.js's output) ----

test('a shipped per-family override (assets.thresholds) is applied without passing explicit params', () => {
  const assets = buildTestAssets({ cog0012: [fixtures.cog0012[0]] }); // only 1 representative, same as above
  assets.thresholds = { defaultParams: {}, familyOverrides: { cog0012: { minRepresentatives: 1 } } };
  const frames = padFrames(fixtures.cog0012[0].seq);
  const calls = searchContigForMarkers(frames, assets); // no explicit params — override must come from assets
  assert.strictEqual(calls.length, 1, 'the shipped override should have let this single-representative call through');
  assert.strictEqual(calls[0].family, 'cog0012');
});

test('a per-family override only affects that family, not others evaluated in the same call', () => {
  const assets = buildTestAssets({ cog0012: [fixtures.cog0012[0]], cog0016: [fixtures.cog0016[0]] });
  assets.thresholds = { defaultParams: {}, familyOverrides: { cog0012: { minRepresentatives: 1 } } };
  // cog0016 has no override, so with only 1 representative each, only cog0012 should be called.
  const framesA = padFrames(fixtures.cog0012[0].seq);
  const callsA = searchContigForMarkers(framesA, assets);
  assert.deepStrictEqual(callsA.map((c) => c.family), ['cog0012']);
});

test('an explicit params argument still wins over a shipped override', () => {
  const assets = buildTestAssets({ cog0012: [fixtures.cog0012[0]] });
  assets.thresholds = { defaultParams: {}, familyOverrides: { cog0012: { minRepresentatives: 1 } } };
  const frames = padFrames(fixtures.cog0012[0].seq);
  // Caller explicitly re-tightens minRepresentatives via familyOverrides in
  // the params argument itself, which searchContigForMarkers spreads last.
  const calls = searchContigForMarkers(frames, assets, { familyOverrides: { cog0012: { minRepresentatives: 2 } } });
  assert.strictEqual(calls.length, 0, 'explicit params should override the shipped threshold asset');
});

report();
