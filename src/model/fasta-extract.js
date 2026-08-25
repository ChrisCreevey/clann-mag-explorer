(function () {
  'use strict';

// Per-MAG FASTA extraction (brief §Export — "built via Blob.slice against
// the first-pass index rather than a second full file read"). Pure
// byte-span arithmetic lives here, testable without a real Blob; the
// actual Blob.slice/arrayBuffer calls are DOM/browser-only glue in
// app.js, same pure-logic/DOM-glue split as fasta-index.js itself.
//
// Known limitation (documented, not engineered around, per
// phase1-investigation.md's Phase 1 framing): a contig whose FASTA lines
// weren't uniformly wrapped (`faiEntry.uniform === false`, e.g. mixed
// line widths) can't be safely re-sliced by a computed byte span — a
// non-uniform contig is skipped from extraction rather than risking a
// corrupt sequence, and callers are told which contigs were skipped so
// they can surface it. A gzip-compressed source's offsets are in the
// *decompressed* logical stream (`faiEntry.sourceCompressed`), so slicing
// needs a decompressed byte buffer, not the raw compressed Blob directly
// — see `planFastaExtraction`'s `needsDecompression` return, which
// callers use to decide whether to decompress first.

/**
 * The exact byte range (start offset + length) covering a contig's
 * sequence lines *as originally wrapped in the source file*, so those
 * raw bytes (newlines and all) can be sliced out and reused verbatim
 * under a new header line — no reformatting needed.
 * @param {{offset:number, length:number, lineBases:number, lineWidth:number, uniform:boolean}} faiEntry
 * @returns {{offset:number, byteLength:number}|null} null if the contig's
 *   wrapping wasn't uniform (can't be safely re-sliced) or has zero length
 */
function computeSequenceByteSpan(faiEntry) {
  if (!faiEntry.uniform || faiEntry.length === 0) return null;
  const { offset, length, lineBases, lineWidth } = faiEntry;
  const lineEndBytes = lineWidth - lineBases;
  const fullLines = Math.floor((length - 1) / lineBases);
  const basesInLastLine = length - fullLines * lineBases;
  const byteLength = fullLines * lineWidth + basesInLastLine + lineEndBytes;
  return { offset, byteLength };
}

/**
 * Builds a per-group extraction plan: for each group (typically a MAG id
 * or bin id), which contigs' byte spans to slice and concatenate, and
 * which contigs had to be skipped.
 * @param {object[]} contigRecords - each needs `id`, `header`, `faiEntry`
 * @param {Map<string, string[]>} contigIdsByGroup - group label -> contig ids
 * @returns {Map<string, {entries:Array<{id:string,header:string,span:{offset:number,byteLength:number}}>,
 *   skippedContigIds:string[]}>} same keys as contigIdsByGroup; entries
 *   preserve the group's given contig order; a group with a source that
 *   needed decompression is unaffected here (offsets stay as-is; it's the
 *   caller's job to slice against a decompressed buffer instead of the
 *   raw Blob when any extracted record has `faiEntry.sourceCompressed`)
 */
function planFastaExtraction(contigRecords, contigIdsByGroup) {
  const recordsById = new Map(contigRecords.map((r) => [r.id, r]));
  const plan = new Map();
  for (const [group, contigIds] of contigIdsByGroup) {
    const entries = [];
    const skippedContigIds = [];
    for (const contigId of contigIds) {
      const record = recordsById.get(contigId);
      if (!record) { skippedContigIds.push(contigId); continue; }
      const span = computeSequenceByteSpan(record.faiEntry);
      if (!span) { skippedContigIds.push(contigId); continue; }
      entries.push({ id: record.id, header: record.header, span });
    }
    plan.set(group, { entries, skippedContigIds });
  }
  return plan;
}

const exportsObj = { computeSequenceByteSpan, planFastaExtraction };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.fastaExtract = exportsObj;
}
})();
