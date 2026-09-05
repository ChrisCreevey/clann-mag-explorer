(function () {
  'use strict';

// On-demand ORF highlighting for the contig evidence panel's Sequence
// section (app.js's initContigSequenceSection): given one contig's raw
// sequence, find substantial open reading frames (stop-to-stop stretches
// above a minimum length, across all six frames) and produce highlighted
// HTML for display. Deliberately not part of the main six-frame-
// translation pass every contig already goes through at load time
// (contig-stats.js/marker-genes.js) — this is a coarser, purely visual
// "does this look like real coding sequence" signal computed fresh for
// one contig only when a student actually clicks it, reusing
// translate.js's frame-translation core rather than re-deriving it.
//
// Nucleotide-coordinate mapping, the part easy to get subtly wrong:
// forward frame f's amino-acid index k covers nucleotides
// [f + 3k, f + 3k + 3) directly. A *reverse* frame's amino-acid index k
// is the reverse complement of the forward window ending 3k+f bases from
// the sequence's end — translateReverseFrameCodes (translate.js) walks
// right-to-left, so increasing k moves *leftward* in forward-sequence
// coordinates. A reverse-frame segment spanning amino-acid indices
// [segStart, segEnd) therefore covers the forward nucleotide range
// [n - f - 3*segEnd, n - f - 3*segStart) — the low/high ends swap
// relative to the forward case precisely because k walks backwards.

const { computeBaseCodes } = (typeof module !== 'undefined' && module.exports)
  ? require('./dna-codes')
  : self.ClannMAG.dnaCodes;
const { translateFrameCodes, translateReverseFrameCodes } = (typeof module !== 'undefined' && module.exports)
  ? require('./translate')
  : self.ClannMAG.translate;

const STOP_CHAR_CODE = 42; // '*'
const X_CHAR_CODE = 88; // 'X' — an ambiguous-base codon, not a stop, but not confidently coding either

const DEFAULT_MIN_AA_LENGTH = 100; // ~300 nt — filters short spurious ORFs random sequence produces plenty of, keeping only substantial, likely-real stretches

/**
 * @param {string} sequence - one contig's raw DNA sequence
 * @param {number} [minAaLength] - shortest amino-acid run counted as an ORF
 * @returns {Array<{start:number, end:number, strand:'fwd'|'rev'}>}
 *   nucleotide ranges (end exclusive) in the *forward*-sequence coordinate
 *   space regardless of strand, so callers never need their own strand-
 *   aware coordinate math
 */
function findOrfRanges(sequence, minAaLength = DEFAULT_MIN_AA_LENGTH) {
  const codes = computeBaseCodes(sequence);
  const n = sequence.length;
  const ranges = [];

  function collectSegments(aaBytes, toNucleotideRange) {
    let segStart = 0;
    for (let k = 0; k <= aaBytes.length; k++) {
      const atBoundary = k === aaBytes.length || aaBytes[k] === STOP_CHAR_CODE || aaBytes[k] === X_CHAR_CODE;
      if (!atBoundary) continue;
      if (k - segStart >= minAaLength) {
        const [start, end] = toNucleotideRange(segStart, k);
        ranges.push({ start, end, strand: toNucleotideRange.strand });
      }
      segStart = k + 1;
    }
  }

  for (let f = 0; f < 3; f++) {
    const toRange = (segStart, segEnd) => [f + 3 * segStart, f + 3 * segEnd];
    toRange.strand = 'fwd';
    collectSegments(translateFrameCodes(codes, f), toRange);
  }
  for (let f = 0; f < 3; f++) {
    const toRange = (segStart, segEnd) => [n - f - 3 * segEnd, n - f - 3 * segStart];
    toRange.strand = 'rev';
    collectSegments(translateReverseFrameCodes(codes, f), toRange);
  }

  return ranges;
}

/**
 * Renders `sequence` as HTML, line-wrapped at `width` characters (`\n`
 * inserted at fixed positions regardless of highlight-span boundaries —
 * valid inside a `<pre>`, since a `<span>` may freely contain a newline),
 * with each ORF range from `findOrfRanges` wrapped in a
 * `.orf-fwd`/`.orf-rev`/`.orf-both` span (styles/main.css) — `.orf-both`
 * only when a forward and reverse ORF genuinely overlap at that position,
 * which real, length-filtered ORFs rarely do but is handled rather than
 * silently favouring one strand. Assumes `sequence` is plain A/C/G/T/N/
 * IUPAC letters (always true for FASTA — no HTML-escaping needed).
 * @param {string} sequence
 * @param {Array<{start:number, end:number, strand:'fwd'|'rev'}>} ranges
 * @param {number} [width]
 * @returns {string} HTML
 */
function highlightOrfsHtml(sequence, ranges, width = 60) {
  const n = sequence.length;
  const marks = new Uint8Array(n); // bit 1 = fwd, bit 2 = rev
  for (const { start, end, strand } of ranges) {
    const bit = strand === 'fwd' ? 1 : 2;
    const from = Math.max(0, start), to = Math.min(n, end);
    for (let i = from; i < to; i++) marks[i] |= bit;
  }
  const classForMark = (m) => (m === 0 ? null : m === 1 ? 'orf-fwd' : m === 2 ? 'orf-rev' : 'orf-both');

  let html = '';
  let i = 0;
  while (i < n) {
    const cls = classForMark(marks[i]);
    let j = i + 1;
    while (j < n && classForMark(marks[j]) === cls) j++;
    let chunk = '';
    for (let p = i; p < j; p++) {
      chunk += sequence[p];
      if ((p + 1) % width === 0 && p + 1 < n) chunk += '\n';
    }
    html += cls ? `<span class="${cls}">${chunk}</span>` : chunk;
    i = j;
  }
  return html;
}

const exportsObj = { DEFAULT_MIN_AA_LENGTH, findOrfRanges, highlightOrfsHtml };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.orfFinder = exportsObj;
}
})();
