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
 * @param {number} [maxSeedsPerSegment] - each stop/X-delimited segment (the
 *   same segmentation this function already gets for free from GROUP_INDEX
 *   rejecting '*'/'X') emits at most roughly this many windows.
 * @param {number} [minStride] - floor on the stride applied to every
 *   segment regardless of length (default 1 = no floor). Real assemblies
 *   are dominated by *short* segments (frequent stop codons in
 *   non-coding-frame noise), so a per-segment cap alone barely touches
 *   their volume — most of the actual lookup-volume reduction has to come
 *   from thinning even short segments uniformly. minStride does that; the
 *   effective stride is max(minStride, ceil(windowCount/maxSeedsPerSegment))
 *   so a segment long enough to exceed the cap at minStride density gets a
 *   further-widened stride to still reach its *full* length rather than
 *   exhausting the budget on its front half (the failure mode a bare fixed
 *   stride has) — this is systematic, evenly-spaced sampling, not random:
 *   for a homologous region of length L anywhere in the segment, evenly
 *   spaced sampling guarantees at least 2 sampled positions fall inside it
 *   once L exceeds ~2x the effective stride, whereas random sampling at the
 *   same density only gets that in expectation and has a worse worst case.
 *   Measured against the real reference index, a query window's index
 *   lookup pulls ~10+ hits on average and every one costs a diagonal-map
 *   get/set — that per-hit bookkeeping, not k-mer generation or extension,
 *   dominates single-threaded runtime, so cutting lookup volume broadly
 *   (minStride) while still guaranteeing full-length coverage on long
 *   segments (the cap/stride interaction) is the highest-leverage
 *   single-threaded perf lever.
 */
function forEachReducedKmer(aaSeq, k, onWindow, maxSeedsPerSegment, minStride) {
  const n = aaSeq.length;
  const cap = maxSeedsPerSegment || Infinity;
  const floorStride = minStride || 1;
  let segStart = -1;
  for (let i = 0; i <= n; i++) {
    const g = i < n ? GROUP_INDEX[aaSeq.charCodeAt(i)] : -1;
    if (g >= 0) {
      if (segStart < 0) segStart = i;
    } else if (segStart >= 0) {
      emitSegmentWindows(aaSeq, segStart, i, k, cap, floorStride, onWindow);
      segStart = -1;
    }
  }
}

/** Emits up to `cap` k-mer windows from aaSeq[start, end), stride at least `floorStride`, spread across the whole segment if that's not enough to fit within `cap`. */
function emitSegmentWindows(aaSeq, start, end, k, cap, floorStride, onWindow) {
  const numWindows = end - start - k + 1;
  if (numWindows <= 0) return;
  const stride = Math.max(floorStride, Math.ceil(numWindows / cap));
  for (let pos = start; pos <= end - k; pos += stride) {
    let code = 0;
    for (let j = 0; j < k; j++) code = code * ALPHABET_SIZE + GROUP_INDEX[aaSeq.charCodeAt(pos + j)];
    onWindow(code, pos);
  }
}

const exportsObj = { MURPHY10_GROUPS, ALPHABET_SIZE, GROUP_INDEX, forEachReducedKmer, emitSegmentWindows };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.reducedAlphabet = exportsObj;
}
})();
