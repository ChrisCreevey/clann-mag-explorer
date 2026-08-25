#!/usr/bin/env node
'use strict';

// Calibrates per-family marker-gene search thresholds against the actual
// shipped index, rather than shipping one global guess for all 40
// families (see docs/phase1-investigation.md "Phase 3 calibration
// findings" for the reasoning and the spot-check that motivated this).
//
// Uses held-out sequences build/01-cluster.js discarded during clustering
// (real, known-family members that never went into the shipped index) as
// free true-positive test data — no external genome download needed. For
// each family: how many representatives does a genuine held-out member of
// that family actually support, and how close does the *next-best*
// competing family's score get? Two families are not equally easy — a
// spot check (COG0012 vs COG0016) found one with ~100x the natural
// representative-agreement headroom of the other, so a single global
// threshold is either too loose for the fragile ones or (if tightened
// globally to be safe) needlessly strict everywhere else.
//
// Only emits an override where the family's own calibration data shows
// the global default would cost it real sensitivity (drop below what its
// own true positives can support) — deliberately not tightening beyond
// the global default even where a family has headroom to spare, since
// that would be fitting to this reference set's *internal* diversity, not
// to real-world divergence (this held-out set is drawn from the same
// collection the index itself came from, so it's a good relative signal
// between families, not a substitute for the brief's own calibration ask
// of testing against genomes not enriched for these families).
//
// Input:  data/scg40-index.bin, data/scg40-refseqs.bin, data/scg40-families.json
//         build/intermediate/scg40_heldout.fasta (from 01-cluster.js)
// Output: data/scg40-thresholds.json
//         build/intermediate/scg40_calibration_report.txt (full per-family
//         diagnostics, for eyeballing what got overridden and why)

const fs = require('fs');
const path = require('path');
const { parseFasta, parseHeader } = require('./02-index');
const { parseAssets, computeFamilyCandidates, DEFAULT_PARAMS } = require('../src/model/marker-genes');

const HELDOUT_INPUT = path.join(__dirname, 'intermediate', 'scg40_heldout.fasta');
const DATA_DIR = path.join(__dirname, '..', 'data');
const REPORT_OUT = path.join(__dirname, 'intermediate', 'scg40_calibration_report.txt');
const THRESHOLDS_OUT = path.join(DATA_DIR, 'scg40-thresholds.json');

const SAMPLE_PER_FAMILY = 150; // evenly-strided sample, not every held-out sequence — build-time cost tradeoff
const SAFETY_PERCENTILE = 0.10; // want ~90% of a family's true positives to clear the derived threshold
const SAFETY_DISCOUNT = 0.8; // shave the percentile value further, since we only sampled, not exhaustively tested

function toArrayBuffer(buf) { return buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length); }

function evenSample(arr, n) {
  if (arr.length <= n) return arr;
  const stride = arr.length / n;
  const out = [];
  for (let i = 0; i < n; i++) out.push(arr[Math.floor(i * stride)]);
  return out;
}

function percentile(sortedAsc, p) {
  if (sortedAsc.length === 0) return null;
  const idx = Math.min(sortedAsc.length - 1, Math.floor(p * sortedAsc.length));
  return sortedAsc[idx];
}

function loadAssets() {
  return parseAssets({
    indexBuf: toArrayBuffer(fs.readFileSync(path.join(DATA_DIR, 'scg40-index.bin'))),
    refSeqsBuf: toArrayBuffer(fs.readFileSync(path.join(DATA_DIR, 'scg40-refseqs.bin'))),
    familyNames: JSON.parse(fs.readFileSync(path.join(DATA_DIR, 'scg40-families.json'), 'utf8')),
  });
}

function main() {
  console.log('Loading shipped assets...');
  const assets = loadAssets();

  console.log(`Reading ${HELDOUT_INPUT}...`);
  const heldOut = parseFasta(fs.readFileSync(HELDOUT_INPUT, 'utf8')).map((r) => ({ ...r, ...parseHeader(r.header) }));
  console.log(`${heldOut.length} held-out sequences.`);

  const byFamily = new Map();
  for (const r of heldOut) {
    if (!byFamily.has(r.family)) byFamily.set(r.family, []);
    byFamily.get(r.family).push(r);
  }

  const reportLines = [
    'family\tsampled\trecall_failures\trepCount_p10\tmargin_p10\trecommended_minReps\trecommended_minMargin',
  ];
  const familyOverrides = {};
  const t0 = Date.now();
  let totalQueries = 0;

  for (const [family, records] of [...byFamily.entries()].sort()) {
    const sample = evenSample(records, SAMPLE_PER_FAMILY);
    const repCounts = [];
    const margins = [];
    let recallFailures = 0;

    for (const r of sample) {
      totalQueries++;
      const frames = [r.seq, '', '', '', '', '']; // held-out sequences are already AA; treat as a single "frame"
      const { byFamily: candidatesByFamily, familyBestScore } = computeFamilyCandidates(frames, assets, DEFAULT_PARAMS);

      const ownHits = candidatesByFamily.get(family);
      if (!ownHits) { recallFailures++; continue; } // didn't even clear the base score/coverage threshold

      const ownBest = familyBestScore.get(family);
      let bestOther = 0;
      for (const [otherFamily, score] of familyBestScore) {
        if (otherFamily !== family && score > bestOther) bestOther = score;
      }
      repCounts.push(ownHits.length);
      margins.push(ownBest - bestOther);
    }

    repCounts.sort((a, b) => a - b);
    margins.sort((a, b) => a - b);
    const repP10 = percentile(repCounts, SAFETY_PERCENTILE);
    const marginP10 = percentile(margins, SAFETY_PERCENTILE);

    const recommendedMinReps = repP10 === null ? null : Math.max(1, Math.floor(repP10 * SAFETY_DISCOUNT));
    const recommendedMinMargin = marginP10 === null ? null : Math.max(0, Math.floor(marginP10 * SAFETY_DISCOUNT));

    const override = {};
    // Only override when the family's own data shows the global default
    // would cost it real sensitivity — never tighten beyond default here
    // (see header note on not overfitting to this reference set's
    // internal diversity).
    if (recommendedMinReps !== null && recommendedMinReps < DEFAULT_PARAMS.minRepresentatives) {
      override.minRepresentatives = recommendedMinReps;
    }
    if (recommendedMinMargin !== null && recommendedMinMargin < DEFAULT_PARAMS.minMargin) {
      override.minMargin = recommendedMinMargin;
    }
    if (Object.keys(override).length > 0) familyOverrides[family] = override;

    reportLines.push(
      `${family}\t${sample.length}\t${recallFailures}\t${repP10 ?? 'n/a'}\t${marginP10 ?? 'n/a'}\t${recommendedMinReps ?? 'n/a'}\t${recommendedMinMargin ?? 'n/a'}`
    );
    console.log(
      `  ${family}: sampled ${sample.length}, ${recallFailures} recall failures, rep p10=${repP10 ?? 'n/a'}, margin p10=${marginP10 ?? 'n/a'}` +
      (override.minRepresentatives || override.minMargin ? ` -> OVERRIDE ${JSON.stringify(override)}` : '')
    );
  }

  fs.writeFileSync(REPORT_OUT, reportLines.join('\n') + '\n');

  const thresholds = { defaultParams: DEFAULT_PARAMS, familyOverrides };
  fs.mkdirSync(DATA_DIR, { recursive: true });
  fs.writeFileSync(THRESHOLDS_OUT, JSON.stringify(thresholds, null, 2));

  console.log(`\nDone in ${((Date.now() - t0) / 1000).toFixed(1)}s, ${totalQueries} total searches.`);
  console.log(`${Object.keys(familyOverrides).length} / ${byFamily.size} families got an override.`);
  console.log(`Wrote ${REPORT_OUT}`);
  console.log(`Wrote ${THRESHOLDS_OUT}`);
}

if (require.main === module) main();

module.exports = { evenSample, percentile };
