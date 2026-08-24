(function () {
  'use strict';

// Per-contig statistics computed on the streaming FASTA pass (brief
// §Per-contig properties / §The FASTA is never held in memory): length,
// GC content, GC skew, a tetranucleotide composition signature, and a
// coding-density estimate — everything retained after the raw sequence
// itself is discarded.

const { translateSixFrames } = (typeof module !== 'undefined' && module.exports)
  ? require('./translate')
  : window.ClannMAG.translate;

const K = 4; // tetranucleotide, the standard choice for binning-style composition vectors
const BASES = ['A', 'C', 'G', 'T'];
const DNA_COMPLEMENT = { A: 'T', C: 'G', G: 'C', T: 'A' };

function reverseComplementDNA(kmer) {
  let out = '';
  for (let i = kmer.length - 1; i >= 0; i--) out += DNA_COMPLEMENT[kmer[i]];
  return out;
}

// Canonical-kmer lookup, built once and cached: every possible 4-mer maps
// to whichever of {itself, its reverse complement} sorts first, collapsing
// the 256 possible 4-mers to 136 canonical bins so strand doesn't matter.
let canonicalKmerCache = null;
function getCanonicalKmerIndex() {
  if (canonicalKmerCache) return canonicalKmerCache;

  const all = [];
  (function gen(prefix) {
    if (prefix.length === K) { all.push(prefix); return; }
    for (const b of BASES) gen(prefix + b);
  })('');

  const canonicalSet = new Set();
  for (const kmer of all) {
    const rc = reverseComplementDNA(kmer);
    canonicalSet.add(kmer < rc ? kmer : rc);
  }
  const canonicalList = [...canonicalSet].sort();
  const canonIndexByKmer = new Map(canonicalList.map((k, i) => [k, i]));

  const kmerToCanonIndex = new Map();
  for (const kmer of all) {
    const rc = reverseComplementDNA(kmer);
    const canon = kmer < rc ? kmer : rc;
    kmerToCanonIndex.set(kmer, canonIndexByKmer.get(canon));
  }

  canonicalKmerCache = { canonicalList, kmerToCanonIndex };
  return canonicalKmerCache;
}

/** Canonical-tetranucleotide relative-frequency signature, keyed by canonical 4-mer. */
function computeTetranucleotideComposition(seq) {
  const { canonicalList, kmerToCanonIndex } = getCanonicalKmerIndex();
  const counts = new Float64Array(canonicalList.length);
  let total = 0;

  for (let i = 0; i + K <= seq.length; i++) {
    const idx = kmerToCanonIndex.get(seq.substr(i, K));
    if (idx === undefined) continue; // window touches a non-ACGT base
    counts[idx]++;
    total++;
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
if (typeof window !== 'undefined') {
  window.ClannMAG = window.ClannMAG || {};
  window.ClannMAG.contigStats = exportsObj;
}
})();
