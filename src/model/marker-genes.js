(function () {
  'use strict';

// Marker-gene identification module (brief §Marker-gene identification
// module). Loads the offline-built seed index, reference sequences, and
// family-name table (build/01-cluster.js -> build/02-index.js — see
// docs/phase1-investigation.md §4-5), and runs the seed-and-extend search
// against a contig's six-frame translation (the same translation
// contig-stats.js already computes for coding density — see
// translate.js's *FromCodes functions).
//
// Scope note: this module produces per-contig family tags + provenance
// taxIDs (the brief's Phase 3 deliverable). The taxID->lineage table and
// LCA/consensus-lineage computation (marker-gene taxonomic consistency,
// brief's per-bin chimerism check) are Phase 6 concerns layered on top of
// these per-contig tags once bins are loaded — build/03-taxonomy.js is
// deferred until that's actually being built, not part of this module.
//
// Calibration placeholders, not validated: X-drop, minScore, minCoverage,
// minMargin, minRepresentatives below are reasonable starting points, not
// numbers tested against negative-control genomes as the brief's Phase 1
// calibration calls for (docs/phase1-investigation.md §6). Exposed as
// `params` so they're easy to override once real calibration happens.

const { forEachReducedKmer } = (typeof module !== 'undefined' && module.exports)
  ? require('./reduced-alphabet')
  : self.ClannMAG.reducedAlphabet;
const { BLOSUM62_SCORE } = (typeof module !== 'undefined' && module.exports)
  ? require('./blosum62')
  : self.ClannMAG.blosum62;

const DEFAULT_PARAMS = {
  xDrop: 15,            // ungapped extension stops once running score falls this far below its best
  minScore: 50,         // minimum total BLOSUM62 alignment score to consider a hit real
  minCoverage: 0.5,      // aligned region must cover at least this fraction of the reference's length
  minMargin: 10,         // best family's score must beat the next-best family's by at least this much
  minRepresentatives: 2,  // at least this many distinct representatives in a family must independently clear threshold
  // Seed-volume perf levers (see reduced-alphabet.js's forEachReducedKmer)
  // — not detection-quality knobs on their own. minSeedStride thins every
  // segment uniformly (real assemblies are dominated by short segments
  // from frequent stop codons in non-coding-frame noise, so this is what
  // actually cuts most of the lookup volume); maxSeedsPerSegment then
  // widens the stride further, but only on segments long enough to
  // otherwise exceed it, so a long ORF still gets sampled across its full
  // length rather than only its front. Swept against real assembly data:
  // minSeedStride 4 + maxSeedsPerSegment 150 was a reasonable cost/recall
  // tradeoff for this tool's "lightweight visualizer, not production"
  // framing (index-lookup/diagonal-map bookkeeping, not extension,
  // dominates single-threaded runtime, so this is the highest-leverage
  // single-threaded lever).
  minSeedStride: 4,
  maxSeedsPerSegment: 150,
  // Gapped-extension prototype (see extendGapped below and
  // docs/scg-blast-verification.md "What's left") — off by default,
  // production behavior is unchanged unless a caller opts in. gapPenalty
  // and bandWidth are unvalidated starting points (a linear gap penalty
  // roughly matching BLOSUM62's scale, and a band wide enough to catch a
  // handful of small indels without either cost exploding or the search
  // wandering too far off the seed diagonal), not calibrated numbers.
  useGappedExtension: false,
  gapPenalty: 8,
  bandWidth: 16,
};

// ---- Binary asset parsing ----

// A dense CSR array indexed directly by k-mer code (size ALPHABET_SIZE^k)
// wastes space proportional to how sparse the reference set's real content
// is — at k=7 only 14.8% of the 10,000,000 possible keys were ever
// populated, so 85% of that array was stored purely to say "no hits here".
// Storing only the populated keys (sorted, CSR-style offsets alongside)
// instead decouples index size from k entirely, at the cost of needing a
// lookup instead of direct array indexing at runtime — buildKeyLookup below
// does that once at parse time (an open-addressed hash table over the
// populated codes), so per-kmer lookup at search time stays O(1) average,
// same as the old direct-index scheme.
function hashCode32(code) {
  return (Math.imul(code ^ (code >>> 16), 0x2545f491) >>> 0);
}

function buildKeyLookup(sortedKeys) {
  const n = sortedKeys.length;
  let capacity = 64;
  while (capacity < n * 2) capacity *= 2; // ~50% load factor
  const mask = capacity - 1;
  const slotIndex = new Int32Array(capacity).fill(-1); // -1 = empty; else index into sortedKeys/keyOffsets
  for (let i = 0; i < n; i++) {
    const code = sortedKeys[i];
    let h = hashCode32(code) & mask;
    while (slotIndex[h] !== -1) h = (h + 1) & mask;
    slotIndex[h] = i;
  }
  return { slotIndex, mask };
}

/** Returns the populated-key index for `code` (for use with keyOffsets), or -1 if code has no hits at all. */
function lookupPopulatedIndex(sortedKeys, keyLookup, code) {
  const { slotIndex, mask } = keyLookup;
  let h = hashCode32(code) & mask;
  for (;;) {
    const idx = slotIndex[h];
    if (idx === -1) return -1;
    if (sortedKeys[idx] === code) return idx;
    h = (h + 1) & mask;
  }
}

// Bloom-filter pre-check for lookupPopulatedIndex: on real assemblies,
// ~97% of query k-mers have no hit in the index at all (see the profiling
// that motivated this — even after every seed-volume cut this session
// made, that's still 30M+ "probably nothing here" lookups per contig
// batch). A Bloom filter never has false negatives, so it's a pure speed
// lever with zero recall risk (unlike every k-mer/index-content change
// earlier in this session) — it just answers "definitely not present"
// cheaper than a full open-addressed probe for the common miss case, at
// the cost of occasionally saying "maybe" for a true non-member (a false
// positive just falls through to the same real lookup as today, so it's
// never wrong, only sometimes not faster). Measured on real assembly data
// before shipping: ~1.1-1.3x on the seeding stage, zero change in results.
//
// Sized via the standard formulas for a target ~1% false-positive rate:
// m ≈ n * 9.585 bits, k ≈ (m/n) * ln2 hash functions — rounded up to a
// power-of-two bit count so lookups can mask instead of modulo.
function buildBloomFilter(sortedKeys) {
  const n = sortedKeys.length;
  let numBits = 64;
  while (numBits < n * 9.585) numBits *= 2;
  const numHashes = Math.max(1, Math.round((numBits / Math.max(n, 1)) * Math.LN2));
  const words = new Uint32Array(Math.ceil(numBits / 32));
  const mask = numBits - 1;
  for (let i = 0; i < n; i++) {
    const code = sortedKeys[i];
    const h1 = hashCode32(code);
    const h2 = (Math.imul(code ^ (code >>> 13), 0x85ebca6b) >>> 0) | 1; // odd step, avoids degenerate stride-0 collisions with h1
    for (let j = 0; j < numHashes; j++) {
      const bit = (h1 + j * h2) >>> 0 & mask;
      words[bit >>> 5] |= (1 << (bit & 31));
    }
  }
  return { words, mask, numHashes };
}

/** true = code might be a populated key (fall through to lookupPopulatedIndex); false = definitely not (skip it). */
function bloomMaybeContains(bloom, code) {
  const { words, mask, numHashes } = bloom;
  const h1 = hashCode32(code);
  const h2 = (Math.imul(code ^ (code >>> 13), 0x85ebca6b) >>> 0) | 1;
  for (let j = 0; j < numHashes; j++) {
    const bit = (h1 + j * h2) >>> 0 & mask;
    if ((words[bit >>> 5] & (1 << (bit & 31))) === 0) return false;
  }
  return true;
}

function parseIndexBinary(buf) {
  const view = new DataView(buf);
  const magic = String.fromCharCode(view.getUint8(0), view.getUint8(1), view.getUint8(2), view.getUint8(3));
  if (magic !== 'SCGI') throw new Error(`parseIndexBinary: bad magic "${magic}", expected "SCGI"`);
  const version = view.getUint8(4);
  if (version !== 2) throw new Error(`parseIndexBinary: unsupported version ${version} (expected 2, the sparse-key format — rebuild via build/02-index.js)`);
  const k = view.getUint8(5);
  const numPopulatedKeys = view.getUint32(8, true);
  const numHits = view.getUint32(12, true);

  let offset = 16;
  const sortedKeys = new Uint32Array(buf, offset, numPopulatedKeys); offset += numPopulatedKeys * 4;
  const keyOffsets = new Uint32Array(buf, offset, numPopulatedKeys + 1); offset += (numPopulatedKeys + 1) * 4;
  const hitRefSeqId = new Uint16Array(buf, offset, numHits); offset += numHits * 2;
  const hitPosition = new Uint16Array(buf, offset, numHits); offset += numHits * 2;

  const keyLookup = buildKeyLookup(sortedKeys);
  const bloom = buildBloomFilter(sortedKeys);

  return { k, numPopulatedKeys, numHits, sortedKeys, keyOffsets, keyLookup, bloom, hitRefSeqId, hitPosition };
}

function parseRefSeqsBinary(buf) {
  const view = new DataView(buf);
  const magic = String.fromCharCode(view.getUint8(0), view.getUint8(1), view.getUint8(2), view.getUint8(3));
  if (magic !== 'SCGR') throw new Error(`parseRefSeqsBinary: bad magic "${magic}", expected "SCGR"`);
  const numSeqs = view.getUint32(8, true);
  const totalResidueBytes = view.getUint32(12, true);

  let offset = 16;
  const seqOffsets = new Uint32Array(buf, offset, numSeqs + 1); offset += (numSeqs + 1) * 4;
  const taxId = new Uint32Array(buf, offset, numSeqs); offset += numSeqs * 4;
  const familyIndex = new Uint8Array(buf, offset, numSeqs); offset += numSeqs;
  const residues = new Uint8Array(buf, offset, totalResidueBytes); offset += totalResidueBytes;

  return { numSeqs, seqOffsets, taxId, familyIndex, residues };
}

/**
 * @param {{indexBuf: ArrayBuffer, refSeqsBuf: ArrayBuffer, familyNames: string[],
 *   thresholds?: {defaultParams?: object, familyOverrides?: Record<string,object>}}} raw
 */
function parseAssets({ indexBuf, refSeqsBuf, familyNames, thresholds }) {
  return { index: parseIndexBinary(indexBuf), refSeqs: parseRefSeqsBinary(refSeqsBuf), familyNames, thresholds };
}

/**
 * Fetch and parse the shipped assets, including the optional per-family
 * calibration table (build/03-calibrate.js's output). Calibration is
 * genuinely optional — the whole marker-gene module already degrades
 * gracefully when assets fail to load (brief's "optional module"
 * framing) — so a missing/failed thresholds fetch just falls back to
 * DEFAULT_PARAMS for every family rather than failing the whole load.
 * @param {string} dataUrl base path, e.g. 'data/'
 */
async function loadMarkerGeneAssets(dataUrl) {
  const [indexBuf, refSeqsBuf, familyNames] = await Promise.all([
    fetch(dataUrl + 'scg40-index.bin').then((r) => r.arrayBuffer()),
    fetch(dataUrl + 'scg40-refseqs.bin').then((r) => r.arrayBuffer()),
    fetch(dataUrl + 'scg40-families.json').then((r) => r.json()),
  ]);
  const thresholds = await fetch(dataUrl + 'scg40-thresholds.json').then((r) => r.json()).catch(() => undefined);
  return parseAssets({ indexBuf, refSeqsBuf, familyNames, thresholds });
}

// ---- Seed-and-extend search ----

/**
 * Ungapped, X-drop extension anchored at (frameAnchor, refAnchor), scored
 * on real residues (BLOSUM62) regardless of what alphabet seeding used.
 * Extends right then left independently from the anchor, each tracking
 * its own running/best score, and sums the two best points — standard
 * BLAST-style ungapped extension. A stop-sentinel or ambiguous residue on
 * either side scores far enough below any real pair (blosum62.js's
 * STOP_SCORE) that extension halts there on its own, with no separate
 * boundary check needed (brief's segmentation-as-side-effect design).
 */
function extendUngapped(frameStr, refResidues, refOffset, refLength, frameAnchor, refAnchor, xDrop) {
  let score = BLOSUM62_SCORE[frameStr.charCodeAt(frameAnchor) * 128 + refResidues[refOffset + refAnchor]];
  let best = score, bestRight = 0;
  for (let i = 1; frameAnchor + i < frameStr.length && refAnchor + i < refLength; i++) {
    score += BLOSUM62_SCORE[frameStr.charCodeAt(frameAnchor + i) * 128 + refResidues[refOffset + refAnchor + i]];
    if (score > best) { best = score; bestRight = i; }
    else if (best - score > xDrop) break;
  }

  let leftScore = 0, bestLeft = 0, bestLeftOffset = 0;
  for (let j = 1; frameAnchor - j >= 0 && refAnchor - j >= 0; j++) {
    leftScore += BLOSUM62_SCORE[frameStr.charCodeAt(frameAnchor - j) * 128 + refResidues[refOffset + refAnchor - j]];
    if (leftScore > bestLeft) { bestLeft = leftScore; bestLeftOffset = j; }
    else if (bestLeft - leftScore > xDrop) break;
  }

  const frameStart = frameAnchor - bestLeftOffset;
  const frameEnd = frameAnchor + bestRight; // inclusive
  return {
    score: best + bestLeft,
    frameStart, frameEnd,
    refStart: refAnchor - bestLeftOffset, refEnd: refAnchor + bestRight,
    alignedLength: frameEnd - frameStart + 1,
  };
}

// ---- Gapped extension (prototype — see docs/scg-blast-verification.md
// "What's left") ----
//
// extendUngapped can't cross an indel: the moment query and reference go
// out of register, every subsequent residue pair looks like noise to
// BLOSUM62, and X-drop halts immediately — regardless of how well the
// alignment resumes a few residues later back in the correct frame. This
// is a banded, X-drop-bounded dynamic-programming extension (the same
// family of algorithm BLAST's own gapped-extension stage uses) that can
// insert a gap to restore register and keep going. Deliberately a LINEAR
// gap penalty for this prototype, not the affine (open+extend) model a
// shipped version would want — simpler to get correct first, and enough
// to test whether crossing indels recovers real recall before investing
// in the fuller model.
//
// Banded: only reference offsets within `bandWidth` of the query offset
// are considered at all, not a full O(query x reference) table — indels
// in real, related sequences are usually small, and this is what keeps
// the extra cost proportional to (extension length x band width) rather
// than to the full reference length. Only ever invoked on (refSeqId,
// diagonal) pairs that already cleared MIN_SEEDS_PER_DIAGONAL, i.e. a
// tiny, pre-filtered slice of all candidates — same cost structure as
// extendUngapped today, just a constant (band-width) factor more
// expensive per call.

/**
 * One-directional banded gapped extension, X-drop bounded. `query`/`ref`
 * are read forward from index 0 (caller reverses them for the leftward
 * direction, matching extendUngapped's left/right split).
 * @returns {{score:number, queryConsumed:number, refConsumed:number}}
 *   the best-scoring prefix pair found, not necessarily consuming
 *   everything available (X-drop/length bounds may cut it short)
 */
function extendGappedOneDirection(query, queryLen, ref, refLen, xDrop, gapPenalty, bandWidth) {
  const width = bandWidth * 2 + 1;
  const NEG_INF = -1e9;
  let prevRow = new Float64Array(width).fill(NEG_INF);
  let currRow = new Float64Array(width).fill(NEG_INF);
  prevRow[bandWidth] = 0; // score[i=0][b=0] = 0 (no residues consumed yet)

  let globalBest = 0, bestQueryConsumed = 0, bestRefConsumed = 0;
  const maxI = queryLen; // i = query residues consumed so far

  for (let i = 1; i <= maxI; i++) {
    currRow.fill(NEG_INF);
    let rowBest = NEG_INF;
    for (let bi = 0; bi < width; bi++) {
      const b = bi - bandWidth; // diagonal offset: j = i + b
      const j = i + b;
      if (j < 0 || j > refLen) continue;

      let cell = NEG_INF;
      if (j >= 1 && prevRow[bi] > NEG_INF) { // diagonal: consume query[i-1] and ref[j-1]
        const sub = BLOSUM62_SCORE[query[i - 1] * 128 + ref[j - 1]];
        const v = prevRow[bi] + sub;
        if (v > cell) cell = v;
      }
      if (bi + 1 < width && prevRow[bi + 1] > NEG_INF) { // gap in reference: consume query[i-1] only (j unchanged)
        const v = prevRow[bi + 1] - gapPenalty;
        if (v > cell) cell = v;
      }
      if (bi - 1 >= 0 && j >= 1 && currRow[bi - 1] > NEG_INF) { // gap in query: consume ref[j-1] only (i unchanged)
        const v = currRow[bi - 1] - gapPenalty;
        if (v > cell) cell = v;
      }

      currRow[bi] = cell;
      if (cell > rowBest) rowBest = cell;
      if (cell > globalBest) { globalBest = cell; bestQueryConsumed = i; bestRefConsumed = j; }
    }
    if (rowBest < globalBest - xDrop) break; // X-drop: every live cell in this row is already a dead end
    const tmp = prevRow; prevRow = currRow; currRow = tmp;
  }

  return { score: globalBest, queryConsumed: bestQueryConsumed, refConsumed: bestRefConsumed };
}

/**
 * Gapped counterpart to extendUngapped, same signature/return shape (so
 * callers can swap between them via params.useGappedExtension without
 * further changes) — `alignedLength` here is the REFERENCE span consumed
 * (frameEnd-frameStart+1 and refEnd-refStart+1 can now differ, since gaps
 * let query/reference spans diverge; callers computing coverage against
 * refLength want the reference span, so that's what's returned as
 * alignedLength for drop-in compatibility with computeFamilyCandidates).
 */
function extendGapped(frameStr, refResidues, refOffset, refLength, frameAnchor, refAnchor, params) {
  const { xDrop, gapPenalty, bandWidth } = params;
  const anchorScore = BLOSUM62_SCORE[frameStr.charCodeAt(frameAnchor) * 128 + refResidues[refOffset + refAnchor]];

  const rightQueryLen = frameStr.length - frameAnchor - 1;
  const rightRefLen = refLength - refAnchor - 1;
  const rightQuery = new Uint8Array(Math.max(0, rightQueryLen));
  for (let i = 0; i < rightQueryLen; i++) rightQuery[i] = frameStr.charCodeAt(frameAnchor + 1 + i);
  const rightRef = refResidues.subarray(refOffset + refAnchor + 1, refOffset + refAnchor + 1 + Math.max(0, rightRefLen));
  const right = extendGappedOneDirection(rightQuery, rightQueryLen, rightRef, rightRefLen, xDrop, gapPenalty, bandWidth);

  const leftQueryLen = frameAnchor;
  const leftRefLen = refAnchor;
  const leftQuery = new Uint8Array(Math.max(0, leftQueryLen));
  for (let i = 0; i < leftQueryLen; i++) leftQuery[i] = frameStr.charCodeAt(frameAnchor - 1 - i); // reversed: index 0 = nearest anchor
  const leftRef = new Uint8Array(Math.max(0, leftRefLen));
  for (let i = 0; i < leftRefLen; i++) leftRef[i] = refResidues[refOffset + refAnchor - 1 - i];
  const left = extendGappedOneDirection(leftQuery, leftQueryLen, leftRef, leftRefLen, xDrop, gapPenalty, bandWidth);

  return {
    score: anchorScore + right.score + left.score,
    frameStart: frameAnchor - left.queryConsumed, frameEnd: frameAnchor + right.queryConsumed,
    refStart: refAnchor - left.refConsumed, refEnd: refAnchor + right.refConsumed,
    alignedLength: left.refConsumed + 1 + right.refConsumed, // reference span, for coverage — see doc comment above
    queryAlignedLength: left.queryConsumed + 1 + right.queryConsumed,
  };
}

// Minimum distinct seeds a (refSeqId, diagonal) must accumulate before
// extension runs — BLAST's classic "two-hit" heuristic. Necessary in
// practice, not just a nicety: measured directly, a single-hit trigger
// made *every* window a candidate extension almost regardless of content
// (81.5% of possible reduced 5-mer keys are populated in this reference
// set — real protein sequences aren't remotely uniform over a 10-letter
// alphabet), so even a fully random, marker-free contig triggered
// extension on the large majority of its windows: ~1.6s for a 20kb random
// contig, which projects to tens of minutes for a realistic assembly.
// A true homologous region produces many seeds on the same diagonal; an
// isolated single seed is usually chance alignment-space noise. Requiring
// 2 collapses random-contig search time by >10x (see
// docs/phase1-investigation.md "Phase 3 findings") while real detection
// is unaffected — a genuine hit clears this by a wide margin.
const MIN_SEEDS_PER_DIAGONAL = 2;

// Flat, reused open-addressed hash table for per-frame diagonal tracking,
// replacing a JS Map<number, object> that used to be allocated fresh per
// frame call. Measured directly: with seeding volume already cut by
// maxSeedsPerSegment/minSeedStride (see reduced-alphabet.js), this
// per-hit bookkeeping — not extension, not k-mer generation — was still
// >60% of single-threaded runtime (231M diagonal lookups across a 3000-
// contig sample), so it's worth a table that avoids both the Map's
// per-entry hashing overhead and per-new-diagonal object allocation/GC.
//
// Sized once, module-level, and "cleared" between frame calls by bumping
// diagEpochCounter rather than by zeroing the arrays — a slot only counts
// as occupied if diagEpoch[slot] matches the current call's epoch, so an
// old call's leftover data is implicitly stale without ever being wiped.
// Capacity (2^20) is a generous multiple of the largest per-frame hit
// count seen on real assemblies, so linear-probing stays cheap (low load
// factor) even for unusually long contigs.
const DIAG_TABLE_CAPACITY = 1 << 20;
const DIAG_TABLE_MASK = DIAG_TABLE_CAPACITY - 1;
const diagTableKeys = new Uint32Array(DIAG_TABLE_CAPACITY);
const diagTableEpoch = new Uint32Array(DIAG_TABLE_CAPACITY);
const diagTableFramePos = new Int32Array(DIAG_TABLE_CAPACITY);
const diagTableRefPos = new Int32Array(DIAG_TABLE_CAPACITY);
const diagTableCount = new Uint8Array(DIAG_TABLE_CAPACITY);
const diagTableExtended = new Uint8Array(DIAG_TABLE_CAPACITY);
let diagTableEpochCounter = 0;

function hashDiagKey(diagKey) {
  return (Math.imul(diagKey ^ (diagKey >>> 15), 0x2545f491) >>> 0) & DIAG_TABLE_MASK;
}

/** Finds diagKey's slot in the current epoch, claiming an empty (stale-epoch) one if it's not already present. */
function findOrCreateDiagSlot(diagKey) {
  let idx = hashDiagKey(diagKey);
  for (let probes = 0; probes < DIAG_TABLE_CAPACITY; probes++) {
    if (diagTableEpoch[idx] !== diagTableEpochCounter) {
      diagTableEpoch[idx] = diagTableEpochCounter;
      diagTableKeys[idx] = diagKey;
      diagTableCount[idx] = 0;
      diagTableExtended[idx] = 0;
      return idx;
    }
    if (diagTableKeys[idx] === diagKey) return idx;
    idx = (idx + 1) & DIAG_TABLE_MASK;
  }
  // Sizing keeps this unreachable in practice (see DIAG_TABLE_CAPACITY note above).
  throw new Error('diagonal table full — unexpectedly high hit volume for one frame');
}

/**
 * Seed one translated frame against the index, extend each qualifying
 * (refSeqId, diagonal) — one that accumulates at least
 * MIN_SEEDS_PER_DIAGONAL distinct seeds — exactly once, and fold results
 * into `bestByRefSeqId` (refSeqId -> best-scoring extension seen across
 * all six frames so far).
 */
function searchFrameAgainstIndex(frameStr, frameIdx, assets, p, bestByRefSeqId) {
  const { k, sortedKeys, keyLookup, bloom, keyOffsets, hitRefSeqId, hitPosition } = assets.index;
  const { seqOffsets, residues } = assets.refSeqs;
  const { xDrop, maxSeedsPerSegment, minSeedStride, useGappedExtension, gapPenalty, bandWidth } = p;
  diagTableEpochCounter = (diagTableEpochCounter + 1) >>> 0;
  if (diagTableEpochCounter === 0) diagTableEpochCounter = 1; // skip the sentinel value on the very unlikely 32-bit wraparound

  forEachReducedKmer(frameStr, k, (code, framePos) => {
    if (!bloomMaybeContains(bloom, code)) return; // definitely no hits — skip the full lookup entirely
    const populatedIdx = lookupPopulatedIndex(sortedKeys, keyLookup, code);
    if (populatedIdx < 0) return; // bloom false positive — code has no hits anywhere in the reference set
    const start = keyOffsets[populatedIdx], end = keyOffsets[populatedIdx + 1];
    for (let h = start; h < end; h++) {
      const refSeqId = hitRefSeqId[h];
      const refPos = hitPosition[h];
      const diagKey = refSeqId * 100000 + (framePos - refPos + 50000);

      const slot = findOrCreateDiagSlot(diagKey);
      if (diagTableCount[slot] === 0) { diagTableFramePos[slot] = framePos; diagTableRefPos[slot] = refPos; }
      if (diagTableCount[slot] < 255) diagTableCount[slot]++; // saturating — only ever compared against MIN_SEEDS_PER_DIAGONAL
      if (diagTableExtended[slot] || diagTableCount[slot] < MIN_SEEDS_PER_DIAGONAL) continue;
      diagTableExtended[slot] = 1;

      const refOffset = seqOffsets[refSeqId];
      const refLength = seqOffsets[refSeqId + 1] - refOffset;
      const result = useGappedExtension
        ? extendGapped(frameStr, residues, refOffset, refLength, diagTableFramePos[slot], diagTableRefPos[slot], { xDrop, gapPenalty, bandWidth })
        : extendUngapped(frameStr, residues, refOffset, refLength, diagTableFramePos[slot], diagTableRefPos[slot], xDrop);

      const prev = bestByRefSeqId.get(refSeqId);
      if (!prev || result.score > prev.score) {
        bestByRefSeqId.set(refSeqId, { ...result, frame: frameIdx, refSeqId });
      }
    }
  }, maxSeedsPerSegment, minSeedStride);
}

/**
 * Run the full search for one contig's six translated frames, applying
 * the three paralog-safety checks together (brief: "not individually,
 * given the explicit risk of mis-assigning paralogs with only positive
 * examples to match against"):
 *  - score/coverage threshold (per hit)
 *  - margin between the best and second-best *family*'s score
 *  - multi-representative agreement within the winning family
 *
 * @param {[string,string,string,string,string,string]} sixFrameTranslations
 * @param {object} assets - from loadMarkerGeneAssets/parseAssets
 * @param {object} [params] - overrides for DEFAULT_PARAMS
 * @returns {Array<{family: string, representativeCount: number, bestScore: number, provenanceTaxId: number}>}
 */
/**
 * Merges DEFAULT_PARAMS, any shipped per-family calibration
 * (assets.thresholds, from build/03-calibrate.js — see
 * docs/phase1-investigation.md "Phase 3 calibration findings"),
 * caller-supplied overrides, then that same family's own override on top
 * — the most specific source wins. Every family without its own override
 * just gets the merged global defaults.
 */
function resolveParams(family, p) {
  const override = p.familyOverrides && p.familyOverrides[family];
  return override ? { ...p, ...override } : p;
}

/**
 * Seeds and extends against every candidate reference sequence, then
 * applies only the per-hit score/coverage threshold (not the margin or
 * multi-representative-agreement checks) — the shared first half of
 * searchContigForMarkers, also used directly by build/03-calibrate.js to
 * inspect the full pre-filter picture (how many representatives a held-out
 * sequence *could* have supported, not just whether it happened to clear
 * whatever the current margin/agreement thresholds are).
 *
 * @returns {{byFamily: Map<string, object[]>, familyBestScore: Map<string, number>}}
 */
function computeFamilyCandidates(sixFrameTranslations, assets, p) {
  const bestByRefSeqId = new Map();
  sixFrameTranslations.forEach((frameStr, frameIdx) => {
    searchFrameAgainstIndex(frameStr, frameIdx, assets, p, bestByRefSeqId);
  });

  const { seqOffsets, familyIndex, taxId } = assets.refSeqs;
  const byFamily = new Map();
  for (const [refSeqId, result] of bestByRefSeqId) {
    const family = assets.familyNames[familyIndex[refSeqId]];
    const fp = resolveParams(family, p);
    const refLength = seqOffsets[refSeqId + 1] - seqOffsets[refSeqId];
    const coverage = result.alignedLength / refLength;
    if (result.score < fp.minScore || coverage < fp.minCoverage) continue;

    if (!byFamily.has(family)) byFamily.set(family, []);
    byFamily.get(family).push({ refSeqId, score: result.score, coverage, family, taxId: taxId[refSeqId] });
  }

  const familyBestScore = new Map();
  for (const [family, hits] of byFamily) familyBestScore.set(family, Math.max(...hits.map((h) => h.score)));

  return { byFamily, familyBestScore };
}

function searchContigForMarkers(sixFrameTranslations, assets, params = {}) {
  const p = {
    ...DEFAULT_PARAMS,
    ...(assets.thresholds && assets.thresholds.defaultParams),
    familyOverrides: (assets.thresholds && assets.thresholds.familyOverrides) || {},
    ...params,
  };
  const { byFamily, familyBestScore } = computeFamilyCandidates(sixFrameTranslations, assets, p);

  const calledFamilies = [];
  for (const [family, hits] of byFamily) {
    const fp = resolveParams(family, p);
    const bestScore = familyBestScore.get(family);
    let secondBest = 0;
    for (const [otherFamily, score] of familyBestScore) {
      if (otherFamily !== family && score > secondBest) secondBest = score;
    }
    if (bestScore - secondBest < fp.minMargin) continue; // fails reciprocal-best-hit-style margin check
    if (hits.length < fp.minRepresentatives) continue; // fails multi-representative agreement

    const best = hits.reduce((a, b) => (b.score > a.score ? b : a));
    calledFamilies.push({
      family,
      representativeCount: hits.length,
      bestScore: best.score,
      provenanceTaxId: best.taxId,
    });
  }
  return calledFamilies;
}

const exportsObj = {
  DEFAULT_PARAMS,
  parseIndexBinary, parseRefSeqsBinary, parseAssets, loadMarkerGeneAssets,
  buildKeyLookup, lookupPopulatedIndex, buildBloomFilter, bloomMaybeContains,
  extendUngapped, extendGapped, extendGappedOneDirection,
  searchFrameAgainstIndex, computeFamilyCandidates, resolveParams, searchContigForMarkers,
};
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.markerGenes = exportsObj;
}
})();
