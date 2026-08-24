#!/usr/bin/env node
'use strict';

// Step 2 of the offline marker-gene build pipeline (see
// docs/phase1-investigation.md §4-5). Builds the seed index from the
// clustered reference set: reduced-alphabet translation (Murphy10, see
// src/model/reduced-alphabet.js — shared with the runtime search so
// windowing is identical on both sides), k=6 windowing, packed into a
// CSR-style binary lookup (data/scg40-index.bin), plus the real
// (non-reduced) reference sequences needed for extension scoring
// (data/scg40-refseqs.bin) and a small family-name sidecar.
//
// k=6, not the brief's lower suggested bound of 5: measured directly (see
// "Phase 3 findings" in docs/phase1-investigation.md), k=5 left 81.5% of
// the 100,000 possible reduced-alphabet keys populated in this reference
// set — real protein sequences aren't close to uniform over a 10-letter
// alphabet, so a random, marker-free contig still seeded on the large
// majority of its windows. k=6 (1,000,000 possible keys) cuts populated-key
// occupancy to ~45%, a real specificity gain, not just a size/speed knob.
//
// Input:  build/intermediate/scg40_clustered.fasta (from 01-cluster.js)
// Output: data/scg40-index.bin
//         data/scg40-refseqs.bin
//         data/scg40-families.json

const fs = require('fs');
const path = require('path');
const { ALPHABET_SIZE, forEachReducedKmer } = require('../src/model/reduced-alphabet');

const K = 6; // reduced-alphabet k-mer length for seeding — see note above on why 6 over the brief's suggested 5
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

function buildIndexBinary(records) {
  const numKeys = ALPHABET_SIZE ** K;
  const keyCounts = new Uint32Array(numKeys);

  // Pass 1: count hits per key so we can compute CSR offsets up front.
  for (let seqId = 0; seqId < records.length; seqId++) {
    forEachReducedKmer(records[seqId].seq, K, (code) => { keyCounts[code]++; });
  }

  const cappedCounts = new Uint32Array(numKeys);
  for (let k = 0; k < numKeys; k++) cappedCounts[k] = Math.min(keyCounts[k], MAX_HITS_PER_KEY);

  const keyOffsets = new Uint32Array(numKeys + 1);
  for (let k = 0; k < numKeys; k++) keyOffsets[k + 1] = keyOffsets[k] + cappedCounts[k];
  const numHits = keyOffsets[numKeys];

  // Uint16 for both: numSeqs and every sequence length are well under
  // 65,536 for this reference set (checked below). Halves the dominant
  // cost in the shipped index (numHits is tens of millions of entries
  // even after capping) versus Uint32.
  if (records.length >= 65536) throw new Error(`buildIndexBinary: ${records.length} ref seqs exceeds Uint16 range — widen hitRefSeqId back to Uint32`);
  const maxLen = Math.max(...records.map((r) => r.seq.length));
  if (maxLen >= 65536) throw new Error(`buildIndexBinary: a ${maxLen}-residue sequence exceeds Uint16 range — widen hitPosition back to Uint32`);

  const hitRefSeqId = new Uint16Array(numHits);
  const hitPosition = new Uint16Array(numHits);

  // Pass 2: gather every hit per key first (needed to sample evenly rather
  // than just keeping the first MAX_HITS_PER_KEY encountered, which would
  // bias toward whichever sequences happen to come first in the file).
  const hitsByKey = new Array(numKeys);
  for (let seqId = 0; seqId < records.length; seqId++) {
    forEachReducedKmer(records[seqId].seq, K, (code, position) => {
      if (!hitsByKey[code]) hitsByKey[code] = [];
      hitsByKey[code].push(seqId, position);
    });
  }

  // Pass 3: place a capped, evenly-strided sample of each key's hits —
  // deterministic (not random) so the shipped index is reproducible
  // across build runs, and spread across the full segment rather than
  // clustered at the start.
  let cursor = 0;
  for (let k = 0; k < numKeys; k++) {
    const hits = hitsByKey[k];
    if (!hits) continue;
    const totalForKey = hits.length / 2;
    const keep = cappedCounts[k];
    const stride = totalForKey / keep;
    for (let i = 0; i < keep; i++) {
      const srcIdx = Math.floor(i * stride);
      hitRefSeqId[cursor] = hits[srcIdx * 2];
      hitPosition[cursor] = hits[srcIdx * 2 + 1];
      cursor++;
    }
  }

  const headerBuf = Buffer.alloc(16);
  headerBuf.write('SCGI', 0, 'ascii');
  headerBuf.writeUInt8(1, 4); // version
  headerBuf.writeUInt8(K, 5);
  headerBuf.writeUInt8(ALPHABET_SIZE, 6);
  headerBuf.writeUInt8(2, 7); // bytes per hitRefSeqId/hitPosition entry (Uint16)
  headerBuf.writeUInt32LE(numKeys, 8);
  headerBuf.writeUInt32LE(numHits, 12);

  return {
    buffer: Buffer.concat([
      headerBuf,
      Buffer.from(keyOffsets.buffer),
      Buffer.from(hitRefSeqId.buffer),
      Buffer.from(hitPosition.buffer),
    ]),
    numHits,
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
  const { buffer: indexBuf, numHits } = buildIndexBinary(records);
  console.log(`  ${numHits.toLocaleString()} hits, ${(indexBuf.length / 1e6).toFixed(1)}MB, ${((Date.now() - t0) / 1000).toFixed(1)}s`);
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
