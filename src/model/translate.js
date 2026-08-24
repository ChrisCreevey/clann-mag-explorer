(function () {
  'use strict';

// Six-frame translation shared by contig-stats.js (coding-density estimate,
// Phase 2) and marker-genes.js (seed-and-extend search, Phase 3) — the
// brief calls for reusing the same translation for both.
//
// Rewritten twice after Phase 2 benchmarking (see docs/phase1-investigation.md
// "Phase 2 findings" / "Performance follow-ups"):
//  1. String concatenation -> typed-array codon lookup.
//  2. Reverse frames -> direct lookup against the forward sequence via
//     RC_CODON_CHAR_CODE, no reverse-complement string ever built.
//  3. (This version) All six frames now read from a single shared
//     Int8Array of 2-bit base codes (dna-codes.js's computeBaseCodes),
//     computed once and also reused by contig-stats.js's GC-counting and
//     composition scan — removing the per-frame charCodeAt+BASE_CODE
//     lookup that six separate string-scanning passes used to repeat.
// Public API (string in, string out) is unchanged for translateFrame/
// translateReverseFrame/translateSixFrames; the *FromCodes variants are
// the new lower-level entry points contig-stats.js calls directly to
// avoid recomputing the code array it already has.

const { BASE_CODE, COMPLEMENT_CHAR_CODE, computeBaseCodes } = (typeof module !== 'undefined' && module.exports)
  ? require('./dna-codes')
  : self.ClannMAG.dnaCodes;

// Standard genetic code (NCBI table 1), keyed by codon string for
// readability, converted once below into a 64-entry numeric lookup
// (index = base2bit0<<4 | base2bit1<<2 | base2bit2). Start-codon
// distinctions (table 11 etc.) don't matter here: nothing downstream
// requires a start codon (see brief's note on contig-edge-truncated
// marker fragments), only the stop/non-stop amino-acid assignment.
const CODON_TABLE = {
  TTT: 'F', TTC: 'F', TTA: 'L', TTG: 'L',
  CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L',
  ATT: 'I', ATC: 'I', ATA: 'I', ATG: 'M',
  GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V',
  TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S',
  CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
  ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
  GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
  TAT: 'Y', TAC: 'Y', TAA: '*', TAG: '*',
  CAT: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q',
  AAT: 'N', AAC: 'N', AAA: 'K', AAG: 'K',
  GAT: 'D', GAC: 'D', GAA: 'E', GAG: 'E',
  TGT: 'C', TGC: 'C', TGA: '*', TGG: 'W',
  CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
  AGT: 'S', AGC: 'S', AGA: 'R', AGG: 'R',
  GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G',
};

const X_CODE = 'X'.charCodeAt(0);

// codonCharCode[b0*16 + b1*4 + b2] = charCode of the amino acid (or '*').
const CODON_CHAR_CODE = new Uint8Array(64);
{
  const code = { A: 0, C: 1, G: 2, T: 3 };
  for (const [codon, aa] of Object.entries(CODON_TABLE)) {
    const idx = (code[codon[0]] << 4) | (code[codon[1]] << 2) | code[codon[2]];
    CODON_CHAR_CODE[idx] = aa.charCodeAt(0);
  }
}

// rcCodonCharCode[b0*16 + b1*4 + b2] = amino acid you'd get by translating
// the *reverse complement* of that forward codon — i.e. the two amino
// acids encoded by a forward triplet and its opposite-strand counterpart
// are both fixed functions of the same 3 bases, so this table makes a
// reverse frame a direct lookup against the forward sequence with no
// reverse-complement string ever built. Reverse complement of codon
// (b0,b1,b2) is (3-b2, 3-b1, 3-b0): reverse the base order, and
// complement each (A/T and C/G pairs sum to 3 under this 2-bit encoding).
const RC_CODON_CHAR_CODE = new Uint8Array(64);
for (let code = 0; code < 64; code++) {
  const b0 = (code >> 4) & 3, b1 = (code >> 2) & 3, b2 = code & 3;
  const rcCode = ((3 - b2) << 4) | ((3 - b1) << 2) | (3 - b0);
  RC_CODON_CHAR_CODE[code] = CODON_CHAR_CODE[rcCode];
}

/** Reverse complement of an uppercase DNA/IUPAC-ambiguity-code string. */
function reverseComplement(seq) {
  const n = seq.length;
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    out[n - 1 - i] = COMPLEMENT_CHAR_CODE[seq.charCodeAt(i)];
  }
  return new TextDecoder().decode(out);
}

/**
 * Translate one reading frame starting at `offset` (0, 1, or 2) against a
 * precomputed base-code array (dna-codes.js's computeBaseCodes), returning
 * raw amino-acid char codes (not yet decoded to a string) — the shared
 * core translateFrame/translateSixFramesFromCodes build on. '*' marks a
 * stop codon (a sentinel, per the brief, that later seeding/extension
 * can't cross — no k-mer lookup key spans one, making ORF segmentation an
 * implicit side effect rather than a separate step); a codon touching any
 * non-ACGT position becomes 'X', distinct from '*' so it never silently
 * acts as a stop boundary.
 */
function translateFrameCodes(codes, offset) {
  const n = Math.max(0, Math.floor((codes.length - offset) / 3));
  const out = new Uint8Array(n);
  for (let k = 0, i = offset; k < n; k++, i += 3) {
    const b0 = codes[i], b1 = codes[i + 1], b2 = codes[i + 2];
    out[k] = (b0 < 0 || b1 < 0 || b2 < 0)
      ? X_CODE
      : CODON_CHAR_CODE[(b0 << 4) | (b1 << 2) | b2];
  }
  return out;
}

/**
 * Reverse-complement reading frame `offset` (same convention as
 * translateFrameCodes, but offset into the reverse complement), computed
 * directly against the forward-oriented `codes` array via
 * RC_CODON_CHAR_CODE — equivalent to
 * translateFrameCodes(reverseComplementCodes(codes), offset) but without
 * ever building that reverse-complement array. Walks `codes` right to
 * left in triplets, matching how the reverse frame's first codon is the
 * reverse complement of the *last* 3-adjusted window of the forward
 * sequence, and each subsequent codon steps left by 3.
 */
function translateReverseFrameCodes(codes, offset) {
  const n = codes.length;
  const count = Math.max(0, Math.floor((n - offset) / 3));
  const out = new Uint8Array(count);
  let end = n - offset; // exclusive end of the current forward-sequence window
  for (let k = 0; k < count; k++) {
    const start = end - 3;
    const b0 = codes[start], b1 = codes[start + 1], b2 = codes[start + 2];
    out[k] = (b0 < 0 || b1 < 0 || b2 < 0)
      ? X_CODE
      : RC_CODON_CHAR_CODE[(b0 << 4) | (b1 << 2) | b2];
    end = start;
  }
  return out;
}

/** String-in/string-out convenience wrapper around translateFrameCodes. */
function translateFrame(seq, offset) {
  return new TextDecoder().decode(translateFrameCodes(computeBaseCodes(seq), offset));
}

/** String-in/string-out convenience wrapper around translateReverseFrameCodes. */
function translateReverseFrame(seq, offset) {
  return new TextDecoder().decode(translateReverseFrameCodes(computeBaseCodes(seq), offset));
}

/**
 * All six reading frames against a precomputed base-code array — the core
 * contig-stats.js calls directly, since it already has the code array
 * from its own GC/composition pass and shouldn't recompute it.
 * @param {Int8Array} codes
 * @returns {[string, string, string, string, string, string]}
 */
function translateSixFramesFromCodes(codes) {
  const decoder = new TextDecoder();
  return [
    decoder.decode(translateFrameCodes(codes, 0)),
    decoder.decode(translateFrameCodes(codes, 1)),
    decoder.decode(translateFrameCodes(codes, 2)),
    decoder.decode(translateReverseFrameCodes(codes, 0)),
    decoder.decode(translateReverseFrameCodes(codes, 1)),
    decoder.decode(translateReverseFrameCodes(codes, 2)),
  ];
}

/**
 * All six reading frames (+1,+2,+3 forward, then -1,-2,-3 on the reverse
 * complement) as one continuous amino-acid string each, per the brief's
 * "no separate ORF-calling step" design.
 *
 * @param {string} seq - DNA sequence (any case)
 * @returns {[string, string, string, string, string, string]}
 */
function translateSixFrames(seq) {
  return translateSixFramesFromCodes(computeBaseCodes(seq));
}

const exportsObj = {
  CODON_TABLE, reverseComplement,
  translateFrame, translateReverseFrame, translateSixFrames,
  translateFrameCodes, translateReverseFrameCodes, translateSixFramesFromCodes,
};
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.translate = exportsObj;
}
})();
