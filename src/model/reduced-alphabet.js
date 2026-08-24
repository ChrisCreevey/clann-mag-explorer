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
 */
function forEachReducedKmer(aaSeq, k, onWindow) {
  let code = 0;
  let validRun = 0;
  const mask = ALPHABET_SIZE ** k;
  for (let i = 0; i < aaSeq.length; i++) {
    const g = GROUP_INDEX[aaSeq.charCodeAt(i)];
    if (g < 0) { validRun = 0; code = 0; continue; }
    code = (code * ALPHABET_SIZE + g) % mask;
    validRun++;
    if (validRun >= k) {
      onWindow(code, i - k + 1);
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
