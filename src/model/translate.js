(function () {
  'use strict';

// Six-frame translation shared by contig-stats.js (coding-density estimate,
// Phase 2) and marker-genes.js (seed-and-extend search, Phase 3) — the
// brief calls for reusing the same translation for both (§Marker-gene
// identification module: "reusing the same six-frame translation, searched
// as continuous sentinel-delimited strings rather than pre-cut segments").
//
// Standard genetic code (NCBI table 1). Start-codon distinctions (table 11
// etc.) don't matter here: nothing downstream requires a start codon (see
// brief's note on contig-edge-truncated marker fragments), only the
// stop/non-stop amino-acid assignment.

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

const DNA_COMPLEMENT = {
  A: 'T', T: 'A', C: 'G', G: 'C', U: 'A',
  R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
  B: 'V', D: 'H', H: 'D', V: 'B', N: 'N',
};

/**
 * Reverse complement of an uppercase DNA/IUPAC-ambiguity-code string.
 * Builds into a pre-sized array and joins once, rather than repeated
 * string `+=`, which matters at assembly scale (tens of thousands of
 * contigs, each translated in six frames — see docs/phase1-investigation.md).
 */
function reverseComplement(seq) {
  const n = seq.length;
  const out = new Array(n);
  for (let i = 0; i < n; i++) {
    out[n - 1 - i] = DNA_COMPLEMENT[seq[i]] || 'N';
  }
  return out.join('');
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
  const out = new Array(n);
  for (let k = 0, i = offset; k < n; k++, i += 3) {
    const codon = seq[i] + seq[i + 1] + seq[i + 2];
    out[k] = CODON_TABLE[codon] || 'X';
  }
  return out.join('');
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
if (typeof window !== 'undefined') {
  window.ClannMAG = window.ClannMAG || {};
  window.ClannMAG.translate = exportsObj;
}
})();
