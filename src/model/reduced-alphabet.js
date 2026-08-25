(function () {
  'use strict';

// Reduced amino-acid alphabet for the marker-gene seed index (brief
// §Marker-gene identification module, step 1: "roughly 10-12 groups of
// biochemically similar residues, collapsing most conservative
// substitutions to the same symbol"). Used both by the offline build
// (build/02-index.js, over reference sequences) and at runtime
// (marker-genes.js, over a contig's six-frame translation) — sharing this
// module keeps the two windowing passes identical, which they must be for
// seeding to find matches at all.
//
// Murphy et al. 2000 ("Simplified amino acid alphabets for protein fold
// recognition...") 10-letter reduced alphabet: LVIM / C / A / G / ST / P /
// FYW / EDNQ / KR / H. A commonly cited grouping, not re-derived here —
// worth revisiting during Phase 1's empirical threshold calibration
// (docs/phase1-investigation.md) if it under/over-collapses in practice.

const MURPHY10_GROUPS = ['LVIM', 'C', 'A', 'G', 'ST', 'P', 'FYW', 'EDNQ', 'KR', 'H'];
const ALPHABET_SIZE = MURPHY10_GROUPS.length; // 10

// charCode -> group index (0-9), -1 for anything not a standard amino
// acid (stop sentinel '*', ambiguous 'X', gaps, whitespace).
const GROUP_INDEX = new Int8Array(128).fill(-1);
MURPHY10_GROUPS.forEach((group, groupIdx) => {
  for (const ch of group) GROUP_INDEX[ch.charCodeAt(0)] = groupIdx;
});

/**
 * Every valid reduced-alphabet k-mer window in `aaSeq`, as a packed
 * base-ALPHABET_SIZE integer code (so it can be used directly as an array
 * index into a CSR-style table). A window touching any non-standard
 * residue (stop sentinel, X, etc.) is skipped — mirrors how the
 * marker-gene search treats a stop-delimited translation frame (brief:
 * "no k-mer lookup key spans a sentinel").
 *
 * @param {string} aaSeq
 * @param {number} k
 * @param {(code: number, position: number) => void} onWindow
 * @param {number} [maxSeedsPerSegment] - once a stop/X-delimited segment
 *   (the same segmentation forEachReducedKmer already gets for free from
 *   GROUP_INDEX rejecting '*'/'X') has emitted this many windows, further
 *   windows in that segment are skipped until the next segment break. A
 *   genuine marker hit anchors near the start of its homologous region, so
 *   this caps the seed volume (and downstream index-lookup/diagonal-map
 *   cost) per long ORF without touching extension, which still runs
 *   against the full segment once a candidate diagonal is found.
 * @param {number} [stride] - only every `stride`-th valid window within a
 *   segment is actually emitted (the skipped windows still update the
 *   rolling hash so the emitted ones are correct, they're just not passed
 *   to onWindow). Measured against the real reference index, a query
 *   window's index lookup pulls ~10+ hits on average and every one costs a
 *   diagonal-map get/set — that per-hit bookkeeping, not k-mer generation
 *   or extension, dominates single-threaded runtime, so cutting the number
 *   of *lookups* (not just per-segment length) is the highest-leverage
 *   single-threaded lever. Safe as long as MIN_SEEDS_PER_DIAGONAL worth of
 *   sampled positions still fall inside any real homologous region — a
 *   region has to be shorter than roughly stride * MIN_SEEDS_PER_DIAGONAL
 *   residues to be missed entirely, well below anything minCoverage would
 *   accept anyway.
 */
function forEachReducedKmer(aaSeq, k, onWindow, maxSeedsPerSegment, stride) {
  let code = 0;
  let validRun = 0;
  let segSeeds = 0;
  let segPos = 0;
  const cap = maxSeedsPerSegment || Infinity;
  const step = stride || 1;
  const mask = ALPHABET_SIZE ** k;
  for (let i = 0; i < aaSeq.length; i++) {
    const g = GROUP_INDEX[aaSeq.charCodeAt(i)];
    if (g < 0) { validRun = 0; code = 0; segSeeds = 0; segPos = 0; continue; }
    code = (code * ALPHABET_SIZE + g) % mask;
    validRun++;
    if (validRun >= k) {
      if (segSeeds < cap && segPos % step === 0) {
        onWindow(code, i - k + 1);
        segSeeds++;
      }
      segPos++;
    }
  }
}

const exportsObj = { MURPHY10_GROUPS, ALPHABET_SIZE, GROUP_INDEX, forEachReducedKmer };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.reducedAlphabet = exportsObj;
}
})();
