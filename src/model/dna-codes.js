(function () {
  'use strict';

// Shared 2-bit base encoding used by translate.js (codon lookup) and
// contig-stats.js (tetranucleotide k-mer scan) to avoid per-base string
// allocation/hashing in their hot loops — see docs/phase1-investigation.md
// "Phase 2 findings" for why this exists (six-frame translation was the
// dominant cost in a 50MB-assembly benchmark).
//
// Loaded via `self` rather than `window` so the same file works unchanged
// on the main thread and inside a Worker (self === window on the main
// thread; Workers have `self` but no `window`).

// charCode -> 2-bit code (A=0,C=1,G=2,T=3), -1 for anything else (N, IUPAC
// ambiguity codes, gaps) so callers can branch on validity without a
// separate character-class check.
const BASE_CODE = new Int8Array(128).fill(-1);
BASE_CODE['A'.charCodeAt(0)] = 0;
BASE_CODE['C'.charCodeAt(0)] = 1;
BASE_CODE['G'.charCodeAt(0)] = 2;
BASE_CODE['T'.charCodeAt(0)] = 3;

const BASE_CHAR = ['A', 'C', 'G', 'T'];

// charCode -> complement charCode, covering standard bases plus IUPAC
// ambiguity codes (unlike BASE_CODE, which only needs the unambiguous 4
// for k-mer/codon indexing). Anything not listed complements to 'N'.
const COMPLEMENT_CHAR_CODE = new Uint8Array(128).fill('N'.charCodeAt(0));
const COMPLEMENT_PAIRS = 'AT TA CG GC UA RY YR SS WW KM MK BV VB DH HD NN';
for (const pair of COMPLEMENT_PAIRS.split(' ')) {
  COMPLEMENT_CHAR_CODE[pair.charCodeAt(0)] = pair.charCodeAt(1);
}

const exportsObj = { BASE_CODE, BASE_CHAR, COMPLEMENT_CHAR_CODE };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.dnaCodes = exportsObj;
}
})();
