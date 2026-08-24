(function () {
  'use strict';

// Six-frame translation shared by contig-stats.js (coding-density estimate,
// Phase 2) and marker-genes.js (seed-and-extend search, Phase 3) — the
// brief calls for reusing the same translation for both.
//
// Rewritten from a string-concatenation implementation to a typed-array
// one after Phase 2 benchmarking showed this was the dominant cost on a
// realistic assembly (see docs/phase1-investigation.md "Phase 2 findings"
// / "Performance follow-ups"): per-codon string allocation
// (`seq[i]+seq[i+1]+seq[i+2]`) and object-key hashing replaced with 2-bit
// base codes, a 64-entry numeric codon table, and a single TextDecoder
// pass instead of per-character string building. Public API (string in,
// string out) is unchanged, so callers and tests didn't need to change.

const { BASE_CODE, COMPLEMENT_CHAR_CODE } = (typeof module !== 'undefined' && module.exports)
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
 * Translate one reading frame starting at `offset` (0, 1, or 2) into a
 * string of one-letter amino acid codes, with '*' marking stop codons —
 * a sentinel, per the brief, that later seeding/extension can't cross
 * (no k-mer lookup key spans a sentinel), making ORF segmentation an
 * implicit side effect rather than a separate step. A codon containing any
 * non-ACGT character (ambiguity codes, gaps) becomes 'X', distinct from
 * '*' so it never silently acts as a stop boundary.
 */
function translateFrame(seq, offset) {
  const n = Math.floor((seq.length - offset) / 3);
  const out = new Uint8Array(n);
  for (let k = 0, i = offset; k < n; k++, i += 3) {
    const b0 = BASE_CODE[seq.charCodeAt(i)];
    const b1 = BASE_CODE[seq.charCodeAt(i + 1)];
    const b2 = BASE_CODE[seq.charCodeAt(i + 2)];
    out[k] = (b0 < 0 || b1 < 0 || b2 < 0)
      ? X_CODE
      : CODON_CHAR_CODE[(b0 << 4) | (b1 << 2) | b2];
  }
  return new TextDecoder().decode(out);
}

/**
 * All six reading frames (+1,+2,+3 forward, then -1,-2,-3 on the reverse
 * complement) as one continuous amino-acid string each, per the brief's
 * "no separate ORF-calling step" design.
 *
 * @param {string} seq - uppercase DNA sequence
 * @returns {[string, string, string, string, string, string]}
 */
function translateSixFrames(seq) {
  const rc = reverseComplement(seq);
  return [
    translateFrame(seq, 0), translateFrame(seq, 1), translateFrame(seq, 2),
    translateFrame(rc, 0), translateFrame(rc, 1), translateFrame(rc, 2),
  ];
}

const exportsObj = { CODON_TABLE, reverseComplement, translateFrame, translateSixFrames };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.translate = exportsObj;
}
})();
