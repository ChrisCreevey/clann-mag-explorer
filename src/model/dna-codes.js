(function () {
  'use strict';

// Shared 2-bit base encoding used by translate.js (codon lookup) and
// contig-stats.js (tetranucleotide k-mer scan + GC counting) to avoid
// per-base string allocation/hashing in their hot loops — see
// docs/phase1-investigation.md "Phase 2 findings" / "Performance
// follow-ups" for the profiling behind this.
//
// Loaded via `self` rather than `window` so the same file works unchanged
// on the main thread and inside a Worker (self === window on the main
// thread; Workers have `self` but no `window`).

// charCode -> 2-bit code (A=0,C=1,G=2,T=3), case-insensitive, -1 for
// anything else (N, IUPAC ambiguity codes, gaps) so callers can branch on
// validity without a separate character-class check. Covering both cases
// here (rather than a separate toUpperCase() pass before lookup) removes
// a full-string-copy pass per contig — see "Performance follow-ups".
const BASE_CODE = new Int8Array(128).fill(-1);
BASE_CODE['A'.charCodeAt(0)] = 0; BASE_CODE['a'.charCodeAt(0)] = 0;
BASE_CODE['C'.charCodeAt(0)] = 1; BASE_CODE['c'.charCodeAt(0)] = 1;
BASE_CODE['G'.charCodeAt(0)] = 2; BASE_CODE['g'.charCodeAt(0)] = 2;
BASE_CODE['T'.charCodeAt(0)] = 3; BASE_CODE['t'.charCodeAt(0)] = 3;

const BASE_CHAR = ['A', 'C', 'G', 'T'];

// charCode -> complement charCode, covering standard bases plus IUPAC
// ambiguity codes (unlike BASE_CODE, which only needs the unambiguous 4
// for k-mer/codon indexing). Anything not listed complements to 'N'.
const COMPLEMENT_CHAR_CODE = new Uint8Array(128).fill('N'.charCodeAt(0));
const COMPLEMENT_PAIRS =
  'AT TA CG GC UA RY YR SS WW KM MK BV VB DH HD NN ' +
  'at ta cg gc ua ry yr ss ww km mk bv vb dh hd nn';
for (const pair of COMPLEMENT_PAIRS.split(' ')) {
  COMPLEMENT_CHAR_CODE[pair.charCodeAt(0)] = pair.charCodeAt(1);
}

/**
 * Convert a raw sequence string into 2-bit base codes, once — the shared
 * starting point for GC counting, tetranucleotide composition, and
 * six-frame translation, which previously each re-derived "which base is
 * this" independently (a separate toUpperCase() copy, string-equality
 * comparisons for GC, and per-frame charCodeAt+lookup for translation).
 * @param {string} seq
 * @returns {Int8Array} same length as seq; -1 marks a non-ACGT position
 */
function computeBaseCodes(seq) {
  const n = seq.length;
  const codes = new Int8Array(n);
  for (let i = 0; i < n; i++) codes[i] = BASE_CODE[seq.charCodeAt(i)];
  return codes;
}

const exportsObj = { BASE_CODE, BASE_CHAR, COMPLEMENT_CHAR_CODE, computeBaseCodes };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.dnaCodes = exportsObj;
}
})();
