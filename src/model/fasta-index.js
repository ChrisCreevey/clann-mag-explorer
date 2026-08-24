(function () {
  'use strict';

// Streaming FASTA reader + .fai-style index builder (brief §The FASTA is
// never held in memory). Reads the given Blob/File exactly once via its
// stream(), one contig's raw sequence buffered at a time (never the whole
// file), and calls `onContig` with a computed per-contig record once each
// record is complete — the raw sequence itself is discarded right after.
//
// Byte-offset tracking assumes single-byte (ASCII) content, true for FASTA
// headers/sequences in normal use; a non-ASCII header would throw the
// offsets off. Detection is content-based (gzip magic bytes), not by file
// name or extension, per the brief's convention elsewhere in this suite.
//
// Known limitation, documented rather than engineered around (brief's own
// framing): for a gzip-compressed input, the resulting .fai-style offsets
// are positions in the *decompressed* logical stream, not in the
// compressed file's bytes — Blob.slice() can't jump straight to them.
// Extraction (Phase 9) against a compressed source needs a full
// re-decompression pass up to that offset; `faiEntry` carries a
// `sourceCompressed` flag so callers can detect this rather than assume
// direct random access always works.

const { computeContigStats } = (typeof module !== 'undefined' && module.exports)
  ? require('./contig-stats')
  : self.ClannMAG.contigStats;

async function isGzipped(blob) {
  if (blob.size < 2) return false;
  const head = new Uint8Array(await blob.slice(0, 2).arrayBuffer());
  return head[0] === 0x1f && head[1] === 0x8b;
}

/**
 * @param {Blob} blob - a File (browser) or Node's web-compatible Blob
 * @param {(record: object) => void} onContig - called once per contig, in
 *   file order, with {id, header, length, gcContent, gcSkew, composition,
 *   codingDensity, ambiguousBaseCount, faiEntry}
 * @returns {Promise<{contigCount: number, totalLength: number, sourceCompressed: boolean}>}
 */
async function streamFasta(blob, onContig) {
  const sourceCompressed = await isGzipped(blob);
  const stream = sourceCompressed
    ? blob.stream().pipeThrough(new DecompressionStream('gzip'))
    : blob.stream();
  const reader = stream.getReader();
  const decoder = new TextDecoder('utf-8');

  let carry = '';
  let byteOffset = 0;
  let current = null;
  let contigCount = 0;
  let totalLength = 0;

  function closeCurrent() {
    if (!current) return;
    const stats = computeContigStats(current.seq);
    const faiEntry = {
      offset: current.firstSeqByteOffset,
      length: stats.length,
      lineBases: current.uniform ? current.lineBases : null,
      lineWidth: current.uniform ? current.lineWidthBytes : null,
      uniform: current.uniform,
      sourceCompressed,
    };
    contigCount++;
    totalLength += stats.length;
    onContig({ id: current.id, header: current.header, ...stats, faiEntry });
    current = null;
  }

  function processLine(rawLine) {
    const hasCR = rawLine.endsWith('\r');
    const line = hasCR ? rawLine.slice(0, -1) : rawLine;
    const lineEndBytes = hasCR ? 2 : 1;

    if (line.startsWith('>')) {
      closeCurrent();
      const header = line.slice(1);
      const id = header.split(/\s+/)[0];
      current = {
        id, header, seq: '',
        firstSeqByteOffset: null, lineBases: null, lineWidthBytes: null,
        uniform: true, sawShortLine: false,
      };
    } else if (current && line.length > 0) {
      current.seq += line;
      if (current.lineBases === null) {
        current.lineBases = line.length;
        current.lineWidthBytes = line.length + lineEndBytes;
        current.firstSeqByteOffset = byteOffset;
      } else if (current.uniform) {
        if (current.sawShortLine || line.length > current.lineBases) {
          current.uniform = false;
        } else if (line.length < current.lineBases) {
          current.sawShortLine = true;
        }
      }
    }
    // Blank lines and any content before the first header are skipped but
    // still counted for byte-offset purposes.
    byteOffset += line.length + lineEndBytes;
  }

  while (true) {
    const { value, done } = await reader.read();
    if (done) break;
    carry += decoder.decode(value, { stream: true });
    const lines = carry.split('\n');
    carry = lines.pop();
    for (const line of lines) processLine(line);
  }
  carry += decoder.decode();
  if (carry.length) processLine(carry);
  closeCurrent();

  return { contigCount, totalLength, sourceCompressed };
}

const exportsObj = { streamFasta, isGzipped };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.fastaIndex = exportsObj;
}
})();
