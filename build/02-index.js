#!/usr/bin/env node
'use strict';

// Step 2 of the offline marker-gene build pipeline (see
// docs/phase1-investigation.md §4-5). Builds the seed index from the
// clustered reference set: reduced-alphabet translation (Murphy10, see
// src/model/reduced-alphabet.js — shared with the runtime search so
// windowing is identical on both sides), k=8 windowing, packed into a
// sparse (populated-keys-only) binary lookup (data/scg40-index.bin), plus
// the real (non-reduced) reference sequences needed for extension scoring
// (data/scg40-refseqs.bin) and a small family-name sidecar.
//
// k history — each step measured directly (see "Phase 3 findings" in
// docs/phase1-investigation.md): k=5 left 81.5% of 100,000 possible
// reduced-alphabet keys populated; k=6, ~45% of 1,000,000; k=7, ~14.8% of
// 10,000,000. Real protein sequences aren't close to uniform over a
// 10-letter alphabet, so lower k means a random, marker-free contig still
// seeds on a large fraction of its windows. k=8 (100,000,000 possible
// keys) continues that trend — ~3% occupancy — for a further real
// specificity gain, not just a size/speed knob: a real query k-mer lands
// on a populated (and, once filled, capped) key far less often, which is
// what actually drives single-threaded runtime (per-hit diagonal
// bookkeeping dominates — see marker-genes.js's diagonal table).
//
// A dense CSR array indexed directly by k-mer code (size ALPHABET_SIZE^K)
// would be unusable at this k — 100,000,000 slots at 4 bytes each is
// 400MB just for the offsets array, before any hit data, regardless of
// how sparse the real content is (and every one of a worker pool's
// threads would hold its own parsed copy). Instead this only stores the
// keys that are actually populated (sorted, with parallel CSR offsets),
// so file size tracks real content rather than 10^K — see
// src/model/marker-genes.js's parseIndexBinary/buildKeyLookup for the
// runtime side (an open-addressed hash table over the populated codes,
// built once at parse time, giving the same O(1)-average lookup a dense
// array would). This also means the *build* itself no longer needs any
// O(ALPHABET_SIZE^K)-sized array (the Map-based counting below), so this
// same code scales to any k without a build-time memory blowup either.
//
// Input:  build/intermediate/scg40_clustered.fasta (from 01-cluster.js)
// Output: data/scg40-index.bin
//         data/scg40-refseqs.bin
//         data/scg40-families.json

const fs = require('fs');
const path = require('path');
const { ALPHABET_SIZE, forEachReducedKmer } = require('../src/model/reduced-alphabet');

const K = 8; // reduced-alphabet k-mer length for seeding — see note above on why 8 over the previous 7
const INPUT = path.join(__dirname, 'intermediate', 'scg40_clustered.fasta');
const DATA_DIR = path.join(__dirname, '..', 'data');

function parseFasta(text) {
  const records = [];
  let header = null, seqParts = null;
  for (const line of text.split('\n')) {
    if (line.startsWith('>')) {
      if (header !== null) records.push({ header, seq: seqParts.join('') });
      header = line.slice(1).trim();
      seqParts = [];
    } else if (seqParts) {
      const trimmed = line.trim();
      if (trimmed) seqParts.push(trimmed);
    }
  }
  if (header !== null) records.push({ header, seq: seqParts.join('') });
  return records;
}

/** Reference header format: COGID.fa.<taxID>.<locus_tag> (brief §Provisioning). */
function parseHeader(header) {
  // Format is COGID.fa.<taxID>.<locus_tag> (brief §Provisioning) — parts[1]
  // is the literal string "fa", the taxID is parts[2].
  const parts = header.split('.');
  const family = parts[0];
  const taxId = Number(parts[2]);
  return { family, taxId };
}

function buildRefSeqsBinary(records, familyIndexByName) {
  const numSeqs = records.length;
  let totalResidueBytes = 0;
  for (const r of records) totalResidueBytes += r.seq.length;

  const seqOffsets = new Uint32Array(numSeqs + 1);
  const taxId = new Uint32Array(numSeqs);
  const familyIndex = new Uint8Array(numSeqs);
  const residues = new Uint8Array(totalResidueBytes);

  let cursor = 0;
  for (let i = 0; i < numSeqs; i++) {
    const r = records[i];
    seqOffsets[i] = cursor;
    for (let j = 0; j < r.seq.length; j++) residues[cursor + j] = r.seq.charCodeAt(j);
    cursor += r.seq.length;
    taxId[i] = r.taxId;
    familyIndex[i] = familyIndexByName.get(r.family);
  }
  seqOffsets[numSeqs] = cursor;

  const headerBuf = Buffer.alloc(16);
  headerBuf.write('SCGR', 0, 'ascii');
  headerBuf.writeUInt8(1, 4); // version
  headerBuf.writeUInt32LE(numSeqs, 8);
  headerBuf.writeUInt32LE(totalResidueBytes, 12);

  return Buffer.concat([
    headerBuf,
    Buffer.from(seqOffsets.buffer),
    Buffer.from(taxId.buffer),
    Buffer.from(familyIndex.buffer),
    Buffer.from(residues.buffer),
  ]);
}

// Cap on hits stored per reduced-alphabet key. The distribution is heavily
// skewed (measured on this reference set at k=6: median 5 hits/key, 90th
// percentile 60) — a handful of generic, low-complexity-ish windows recur
// across many sequences and dominate storage while contributing little
// discriminating power on their own (real specificity comes from several
// hits landing on a consistent diagonal — see marker-genes.js's two-hit
// requirement). Same idea as BLAST-family "neighborhood word" hit-list
// capping. 20 keeps every key at/below the cap (the large majority —
// median is 5) completely untouched, while directly cutting the
// runtime cost of seeding against exactly the generic/over-represented
// keys that would otherwise dominate a random, marker-free contig's
// search time (see "Phase 3 findings" in docs/phase1-investigation.md —
// this cap plus k=6 plus the two-hit heuristic together cut a 50kb random
// contig's search time from ~1.6s to ~0.24s).
const MAX_HITS_PER_KEY = 20;

// Floor on hits stored per reduced-alphabet key — a key seen only this many
// times or fewer across the *entire* 39,854-sequence reference set is
// dropped from the index entirely. Measured at k=8: 57% of populated keys
// are singletons (exactly 1 raw occurrence anywhere in the reference set),
// yet a singleton key can only ever seed *one* reference sequence at *one*
// diagonal — it can never by itself satisfy marker-genes.js's two-hit-per-
// diagonal requirement, so it only ever matters when a second, independent
// seed (very plausibly from a different, non-singleton key nearby) also
// lands on the same diagonal. Dropping singletons (MIN_HITS_PER_KEY=2) cuts
// the populated-key count roughly in half — a real reduction in shipped
// index size — for a recall cost that should be small precisely because
// singleton keys are the ones least likely to be load-bearing on their own.
// Checked against real assembly data (not just reasoned about) before
// shipping — see the commit message for the actual measured recall delta.
const MIN_HITS_PER_KEY = 2;

function buildIndexBinary(records) {
  // Pass 1: gather every hit per key, keyed by a plain Map rather than a
  // dense ALPHABET_SIZE^K-sized array — memory tracks actual populated-key
  // count (millions, not hundreds of millions), so this scales to any k.
  // Collecting every hit (not just the first MAX_HITS_PER_KEY seen) is
  // needed to sample evenly rather than biasing toward whichever sequences
  // happen to come first in the file.
  const hitsByKey = new Map();
  for (let seqId = 0; seqId < records.length; seqId++) {
    forEachReducedKmer(records[seqId].seq, K, (code, position) => {
      let hits = hitsByKey.get(code);
      if (!hits) { hits = []; hitsByKey.set(code, hits); }
      hits.push(seqId, position);
    });
  }

  const sortedCodes = [...hitsByKey.keys()]
    .filter((code) => hitsByKey.get(code).length / 2 >= MIN_HITS_PER_KEY)
    .sort((a, b) => a - b);
  const numPopulatedKeys = sortedCodes.length;
  const sortedKeys = Uint32Array.from(sortedCodes);

  const cappedCounts = new Uint32Array(numPopulatedKeys);
  for (let i = 0; i < numPopulatedKeys; i++) cappedCounts[i] = Math.min(hitsByKey.get(sortedCodes[i]).length / 2, MAX_HITS_PER_KEY);

  const keyOffsets = new Uint32Array(numPopulatedKeys + 1);
  for (let i = 0; i < numPopulatedKeys; i++) keyOffsets[i + 1] = keyOffsets[i] + cappedCounts[i];
  const numHits = keyOffsets[numPopulatedKeys];

  // Uint16 for both: numSeqs and every sequence length are well under
  // 65,536 for this reference set (checked below). Halves the dominant
  // cost in the shipped index (numHits is millions of entries even after
  // capping) versus Uint32.
  if (records.length >= 65536) throw new Error(`buildIndexBinary: ${records.length} ref seqs exceeds Uint16 range — widen hitRefSeqId back to Uint32`);
  const maxLen = Math.max(...records.map((r) => r.seq.length));
  if (maxLen >= 65536) throw new Error(`buildIndexBinary: a ${maxLen}-residue sequence exceeds Uint16 range — widen hitPosition back to Uint32`);

  const hitRefSeqId = new Uint16Array(numHits);
  const hitPosition = new Uint16Array(numHits);

  // Pass 2: place a capped, evenly-strided sample of each key's hits —
  // deterministic (not random) so the shipped index is reproducible
  // across build runs, and spread across the full segment rather than
  // clustered at the start.
  let cursor = 0;
  for (let i = 0; i < numPopulatedKeys; i++) {
    const hits = hitsByKey.get(sortedCodes[i]);
    const totalForKey = hits.length / 2;
    const keep = cappedCounts[i];
    const stride = totalForKey / keep;
    for (let j = 0; j < keep; j++) {
      const srcIdx = Math.floor(j * stride);
      hitRefSeqId[cursor] = hits[srcIdx * 2];
      hitPosition[cursor] = hits[srcIdx * 2 + 1];
      cursor++;
    }
  }

  const headerBuf = Buffer.alloc(16);
  headerBuf.write('SCGI', 0, 'ascii');
  headerBuf.writeUInt8(2, 4); // version 2: sparse (populated-keys-only) format
  headerBuf.writeUInt8(K, 5);
  headerBuf.writeUInt8(ALPHABET_SIZE, 6);
  headerBuf.writeUInt8(2, 7); // bytes per hitRefSeqId/hitPosition entry (Uint16)
  headerBuf.writeUInt32LE(numPopulatedKeys, 8);
  headerBuf.writeUInt32LE(numHits, 12);

  return {
    buffer: Buffer.concat([
      headerBuf,
      Buffer.from(sortedKeys.buffer),
      Buffer.from(keyOffsets.buffer),
      Buffer.from(hitRefSeqId.buffer),
      Buffer.from(hitPosition.buffer),
    ]),
    numHits,
    numPopulatedKeys,
  };
}

function main() {
  console.log(`Reading ${INPUT}...`);
  const records = parseFasta(fs.readFileSync(INPUT, 'utf8')).map((r) => ({ ...r, ...parseHeader(r.header) }));
  console.log(`${records.length} clustered representatives.`);

  const familyNames = [...new Set(records.map((r) => r.family))].sort();
  const familyIndexByName = new Map(familyNames.map((f, i) => [f, i]));
  console.log(`${familyNames.length} families.`);

  fs.mkdirSync(DATA_DIR, { recursive: true });

  console.log('Building seed index...');
  const t0 = Date.now();
  const { buffer: indexBuf, numHits, numPopulatedKeys } = buildIndexBinary(records);
  console.log(`  ${numPopulatedKeys.toLocaleString()} populated keys (of ${(ALPHABET_SIZE ** K).toLocaleString()} possible), ${numHits.toLocaleString()} hits, ${(indexBuf.length / 1e6).toFixed(1)}MB, ${((Date.now() - t0) / 1000).toFixed(1)}s`);
  fs.writeFileSync(path.join(DATA_DIR, 'scg40-index.bin'), indexBuf);

  console.log('Building reference sequence store...');
  const refSeqsBuf = buildRefSeqsBinary(records, familyIndexByName);
  console.log(`  ${(refSeqsBuf.length / 1e6).toFixed(1)}MB`);
  fs.writeFileSync(path.join(DATA_DIR, 'scg40-refseqs.bin'), refSeqsBuf);

  fs.writeFileSync(path.join(DATA_DIR, 'scg40-families.json'), JSON.stringify(familyNames));

  console.log('Done.');
  console.log(`Wrote ${path.join(DATA_DIR, 'scg40-index.bin')}`);
  console.log(`Wrote ${path.join(DATA_DIR, 'scg40-refseqs.bin')}`);
  console.log(`Wrote ${path.join(DATA_DIR, 'scg40-families.json')}`);
}

if (require.main === module) main();

module.exports = { parseFasta, parseHeader, buildIndexBinary, buildRefSeqsBinary };
