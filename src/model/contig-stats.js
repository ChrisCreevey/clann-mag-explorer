(function () {
  'use strict';

// Per-contig statistics computed on the streaming FASTA pass (brief
// §Per-contig properties / §The FASTA is never held in memory): length,
// GC content, GC skew, a tetranucleotide composition signature, and a
// coding-density estimate — everything retained after the raw sequence
// itself is discarded.
//
// The composition scan was rewritten from per-position substr()+Map
// lookup to a rolling 2-bit-per-base integer code after Phase 2
// benchmarking flagged it as a secondary cost alongside six-frame
// translation (see docs/phase1-investigation.md "Performance follow-ups").

const { translateSixFrames } = (typeof module !== 'undefined' && module.exports)
  ? require('./translate')
  : self.ClannMAG.translate;
const { BASE_CODE, BASE_CHAR } = (typeof module !== 'undefined' && module.exports)
  ? require('./dna-codes')
  : self.ClannMAG.dnaCodes;

const K = 4; // tetranucleotide, the standard choice for binning-style composition vectors
const KMER_SPACE = 1 << (2 * K); // 256 possible 4-mers

function codeToKmerString(code) {
  let s = '';
  for (let shift = 2 * (K - 1); shift >= 0; shift -= 2) s += BASE_CHAR[(code >> shift) & 3];
  return s;
}

/** Reverse complement of a k-mer packed as a 2-bit-per-base integer. */
function reverseComplementCode(code, k) {
  let rc = 0;
  for (let i = 0; i < k; i++) {
    const base = (code >> (2 * i)) & 3; // LSB-first walk = kmer read back-to-front
    rc = (rc << 2) | (3 - base); // A(0)<->T(3), C(1)<->G(2): complement is 3-x
  }
  return rc;
}

// Canonical-kmer lookup, built once and cached: every possible 4-mer code
// maps to the array index of whichever of {itself, its reverse complement}
// is numerically smaller, collapsing the 256 possible 4-mers to 136
// canonical bins so strand doesn't matter.
let canonicalKmerCache = null;
function getCanonicalKmerIndex() {
  if (canonicalKmerCache) return canonicalKmerCache;

  const canonicalOfCode = new Int16Array(KMER_SPACE);
  for (let code = 0; code < KMER_SPACE; code++) {
    const rc = reverseComplementCode(code, K);
    canonicalOfCode[code] = Math.min(code, rc);
  }
  const canonicalCodes = [...new Set(canonicalOfCode)].sort((a, b) => a - b);
  const canonIndexByCode = new Map(canonicalCodes.map((c, i) => [c, i]));

  const kmerCodeToCanonIndex = new Int16Array(KMER_SPACE);
  for (let code = 0; code < KMER_SPACE; code++) {
    kmerCodeToCanonIndex[code] = canonIndexByCode.get(canonicalOfCode[code]);
  }
  const canonicalList = canonicalCodes.map(codeToKmerString);

  canonicalKmerCache = { canonicalList, kmerCodeToCanonIndex };
  return canonicalKmerCache;
}

/** Canonical-tetranucleotide relative-frequency signature, keyed by canonical 4-mer. */
function computeTetranucleotideComposition(seq) {
  const { canonicalList, kmerCodeToCanonIndex } = getCanonicalKmerIndex();
  const counts = new Float64Array(canonicalList.length);
  let total = 0;

  let code = 0;
  let validRun = 0; // consecutive ACGT bases seen so far, resets on any ambiguity code
  for (let i = 0; i < seq.length; i++) {
    const b = BASE_CODE[seq.charCodeAt(i)];
    if (b < 0) { validRun = 0; continue; }
    code = ((code << 2) | b) & (KMER_SPACE - 1);
    validRun++;
    if (validRun >= K) {
      counts[kmerCodeToCanonIndex[code]]++;
      total++;
    }
  }

  const freq = {};
  for (let i = 0; i < canonicalList.length; i++) {
    freq[canonicalList[i]] = total ? counts[i] / total : 0;
  }
  return freq;
}

/**
 * Coding-density estimate: fraction of six-frame-translated amino acid
 * positions that fall in a stop-to-stop segment at least `minSegmentLen`
 * residues long — a proxy for "how much of this contig looks like coding
 * sequence", not a real gene call. Segments below the threshold (frequent
 * in non-coding/intergenic stretches, which hit stop codons roughly every
 * ~20 codons by chance alone) don't count as coding.
 */
function computeCodingDensity(frames, minSegmentLen = 20) {
  let codingAaCount = 0;
  let totalAaCount = 0;
  for (const frame of frames) {
    for (const segment of frame.split('*')) {
      totalAaCount += segment.length; // denominator excludes stop-codon sentinels themselves
      if (segment.length >= minSegmentLen) codingAaCount += segment.length;
    }
  }
  return totalAaCount ? codingAaCount / totalAaCount : 0;
}

/**
 * @param {string} seq - raw contig sequence (any case, may contain
 *   ambiguity codes) for this pass only — callers discard it right after.
 */
function computeContigStats(seq) {
  const upper = seq.toUpperCase();
  const length = upper.length;

  let g = 0, c = 0, otherCount = 0;
  for (let i = 0; i < length; i++) {
    const ch = upper[i];
    if (ch === 'G') g++;
    else if (ch === 'C') c++;
    else if (ch !== 'A' && ch !== 'T') otherCount++;
  }
  const gcContent = length ? (g + c) / length : 0;
  const gcSkew = (g + c) ? (g - c) / (g + c) : 0;

  const composition = computeTetranucleotideComposition(upper);
  const frames = translateSixFrames(upper);
  const codingDensity = computeCodingDensity(frames);

  return { length, gcContent, gcSkew, composition, codingDensity, ambiguousBaseCount: otherCount };
}

const exportsObj = { computeContigStats, computeTetranucleotideComposition, computeCodingDensity, getCanonicalKmerIndex };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.contigStats = exportsObj;
}
})();
