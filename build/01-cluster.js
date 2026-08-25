#!/usr/bin/env node
'use strict';

// Step 1 of the offline marker-gene build pipeline (see
// docs/phase1-investigation.md §4). Pre-clusters reference-data/scg40_raw.fasta
// within each COG family at ~90% identity via a fast approximate method:
// per-sequence 5-mer MinHash signatures + greedy single-pass clustering
// (longest sequence first, first-matching-cluster wins — same convention
// as cd-hit), rather than exact all-pairs alignment across ~1,700
// sequences/family x 40 families, which is unnecessary for a 90% cutoff
// and would be far slower for no real benefit at this threshold.
//
// MinHash Jaccard estimate is a proxy for identity, not identity itself —
// flagged in docs/phase1-investigation.md as an open point ("How
// aggressively to pre-cluster...") worth revisiting once tested against
// real student assemblies. The threshold below is a starting point.
//
// Input:  reference-data/scg40_raw.fasta
// Output: build/intermediate/scg40_clustered.fasta (one representative
//         sequence per cluster, header unchanged)
//         build/intermediate/scg40_heldout.fasta (every sequence that
//         clustered INTO an existing representative rather than becoming
//         one itself — i.e. never shipped in the index — used by
//         build/03-calibrate.js as free, known-family held-out test data:
//         a real member of the family the search never got to see)
//         build/intermediate/scg40_cluster_report.txt (per-family
//         original count -> cluster count, for eyeballing how much
//         diversity survives)

const fs = require('fs');
const path = require('path');

const INPUT = path.join(__dirname, '..', 'reference-data', 'scg40_raw.fasta');
const OUT_DIR = path.join(__dirname, 'intermediate');
const OUT_FASTA = path.join(OUT_DIR, 'scg40_clustered.fasta');
const OUT_HELDOUT_FASTA = path.join(OUT_DIR, 'scg40_heldout.fasta');
const OUT_REPORT = path.join(OUT_DIR, 'scg40_cluster_report.txt');

const K = 5; // amino acid k-mer length for the MinHash signature
const NUM_HASHES = 24; // signature size: more = more accurate Jaccard estimate, slower
// MinHash Jaccard estimate proxy for ~90% identity. Rough justification:
// for two same-length sequences differing at a fraction `m` of positions
// (independent substitutions), each mismatch can destroy up to K
// overlapping k-mers, so the expected shared-k-mer fraction has a floor
// around `1 - K*m`. At 90% identity (m=0.1, K=5) that's ~0.5 — matched
// empirically against this reference set (0.5 keeps ~61% of sequences per
// family, 0.3 keeps ~40%, 0.7 keeps ~79%; see build/intermediate/scg40_cluster_report.txt
// after running). Still a proxy, not measured identity — an open point
// for Phase 1 recalibration per docs/phase1-investigation.md.
const SIMILARITY_THRESHOLD = 0.5;

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

// Deterministic 32-bit string hash (FNV-1a), used both to turn a k-mer
// into a number and, seeded per hash-function index, as the basis for
// each MinHash permutation.
function fnv1a(str, seed) {
  let h = 0x811c9dc5 ^ seed;
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  return h >>> 0;
}

function kmerSet(seq) {
  const set = new Set();
  for (let i = 0; i + K <= seq.length; i++) set.add(seq.slice(i, i + K));
  return set;
}

/** MinHash signature: NUM_HASHES independent min-hash values over the k-mer set. */
function minHashSignature(kmers) {
  const sig = new Uint32Array(NUM_HASHES).fill(0xffffffff);
  for (const kmer of kmers) {
    for (let h = 0; h < NUM_HASHES; h++) {
      const v = fnv1a(kmer, h * 2654435761);
      if (v < sig[h]) sig[h] = v;
    }
  }
  return sig;
}

function estimateJaccard(sigA, sigB) {
  let matches = 0;
  for (let i = 0; i < NUM_HASHES; i++) if (sigA[i] === sigB[i]) matches++;
  return matches / NUM_HASHES;
}

function clusterFamily(records) {
  // Longest first: a longer sequence makes a more informative cluster
  // representative (same convention cd-hit uses), and short fragments
  // are more likely to nest inside a longer relative's cluster.
  const ordered = [...records].sort((a, b) => b.seq.length - a.seq.length);
  const representatives = []; // [{ record, signature }]
  const heldOut = [];

  for (const record of ordered) {
    const sig = minHashSignature(kmerSet(record.seq));
    let placed = false;
    for (const rep of representatives) {
      if (estimateJaccard(sig, rep.signature) >= SIMILARITY_THRESHOLD) {
        placed = true;
        break;
      }
    }
    if (placed) heldOut.push(record);
    else representatives.push({ record, signature: sig });
  }
  return { representatives: representatives.map((r) => r.record), heldOut };
}

function main() {
  console.log(`Reading ${INPUT}...`);
  const text = fs.readFileSync(INPUT, 'utf8');
  const records = parseFasta(text);
  console.log(`Parsed ${records.length} sequences.`);

  const byFamily = new Map();
  for (const r of records) {
    const family = r.header.split('.')[0];
    if (!byFamily.has(family)) byFamily.set(family, []);
    byFamily.get(family).push(r);
  }
  console.log(`${byFamily.size} families.`);

  fs.mkdirSync(OUT_DIR, { recursive: true });
  const reportLines = ['family\toriginal_count\tcluster_count\tretained_fraction'];
  const outParts = [];
  const heldOutParts = [];
  let totalOriginal = 0, totalClusters = 0, totalHeldOut = 0;
  const t0 = Date.now();

  for (const [family, famRecords] of [...byFamily.entries()].sort()) {
    const { representatives: reps, heldOut } = clusterFamily(famRecords);
    totalOriginal += famRecords.length;
    totalClusters += reps.length;
    totalHeldOut += heldOut.length;
    reportLines.push(
      `${family}\t${famRecords.length}\t${reps.length}\t${(reps.length / famRecords.length).toFixed(3)}`
    );
    for (const r of reps) outParts.push(`>${r.header}\n${r.seq}\n`);
    for (const r of heldOut) heldOutParts.push(`>${r.header}\n${r.seq}\n`);
    console.log(`  ${family}: ${famRecords.length} -> ${reps.length} (${((reps.length / famRecords.length) * 100).toFixed(1)}%)`);
  }

  fs.writeFileSync(OUT_FASTA, outParts.join(''));
  fs.writeFileSync(OUT_HELDOUT_FASTA, heldOutParts.join(''));
  reportLines.push(`TOTAL\t${totalOriginal}\t${totalClusters}\t${(totalClusters / totalOriginal).toFixed(3)}`);
  fs.writeFileSync(OUT_REPORT, reportLines.join('\n') + '\n');

  console.log(`Done in ${((Date.now() - t0) / 1000).toFixed(1)}s.`);
  console.log(`${totalOriginal} -> ${totalClusters} representatives (${((totalClusters / totalOriginal) * 100).toFixed(1)}%), ${totalHeldOut} held out.`);
  console.log(`Wrote ${OUT_FASTA}`);
  console.log(`Wrote ${OUT_HELDOUT_FASTA}`);
  console.log(`Wrote ${OUT_REPORT}`);
}

main();
