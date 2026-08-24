(function () {
  'use strict';

// Streaming FASTA reader + .fai-style index builder. See
// docs/phase1-investigation.md §3 for the plan this implements. The
// original File/Blob is never fully read into memory: this module reads it
// once, chunk by chunk, retaining only the per-contig numeric record and
// index entry (byte offset, sequence length, line-layout) needed for later
// random-access extraction via blob.slice().
//
// Stubbed pending Phase 2 (per suggested build phases in the brief) —
// implemented after Phase 1's investigation doc is confirmed against a
// real multi-hundred-MB assembly.

/**
 * @param {File} file
 * @param {(record: object) => void} onContig - called once per contig with
 *   {id, length, gc, gcSkew, faiEntry, ...} once its stats are computed and
 *   its sequence discarded.
 */
async function streamFasta(file, onContig) {
  throw new Error('streamFasta: not implemented yet');
}

const exportsObj = { streamFasta };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof window !== 'undefined') {
  window.ClannMAG = window.ClannMAG || {};
  window.ClannMAG.fastaIndex = exportsObj;
}
})();
