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
};

// ---- Binary asset parsing ----

function parseIndexBinary(buf) {
  const view = new DataView(buf);
  const magic = String.fromCharCode(view.getUint8(0), view.getUint8(1), view.getUint8(2), view.getUint8(3));
  if (magic !== 'SCGI') throw new Error(`parseIndexBinary: bad magic "${magic}", expected "SCGI"`);
  const k = view.getUint8(5);
  const numKeys = view.getUint32(8, true);
  const numHits = view.getUint32(12, true);

  let offset = 16;
  const keyOffsets = new Uint32Array(buf, offset, numKeys + 1); offset += (numKeys + 1) * 4;
  const hitRefSeqId = new Uint16Array(buf, offset, numHits); offset += numHits * 2;
  const hitPosition = new Uint16Array(buf, offset, numHits); offset += numHits * 2;

  return { k, numKeys, numHits, keyOffsets, hitRefSeqId, hitPosition };
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

/**
 * Seed one translated frame against the index, extend each qualifying
 * (refSeqId, diagonal) — one that accumulates at least
 * MIN_SEEDS_PER_DIAGONAL distinct seeds — exactly once, and fold results
 * into `bestByRefSeqId` (refSeqId -> best-scoring extension seen across
 * all six frames so far).
 */
function searchFrameAgainstIndex(frameStr, frameIdx, assets, xDrop, bestByRefSeqId, maxSeedsPerSegment, minSeedStride) {
  const { k, keyOffsets, hitRefSeqId, hitPosition } = assets.index;
  const { seqOffsets, residues } = assets.refSeqs;
  const diagonals = new Map(); // diagKey -> { count, framePos, refPos, extended }

  forEachReducedKmer(frameStr, k, (code, framePos) => {
    const start = keyOffsets[code], end = keyOffsets[code + 1];
    for (let h = start; h < end; h++) {
      const refSeqId = hitRefSeqId[h];
      const refPos = hitPosition[h];
      const diagKey = refSeqId * 100000 + (framePos - refPos + 50000);

      let diag = diagonals.get(diagKey);
      if (!diag) {
        diag = { count: 0, framePos, refPos, extended: false };
        diagonals.set(diagKey, diag);
      }
      diag.count++;
      if (diag.extended || diag.count < MIN_SEEDS_PER_DIAGONAL) continue;
      diag.extended = true;

      const refOffset = seqOffsets[refSeqId];
      const refLength = seqOffsets[refSeqId + 1] - refOffset;
      const result = extendUngapped(frameStr, residues, refOffset, refLength, diag.framePos, diag.refPos, xDrop);

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
    searchFrameAgainstIndex(frameStr, frameIdx, assets, p.xDrop, bestByRefSeqId, p.maxSeedsPerSegment, p.minSeedStride);
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
  extendUngapped, searchFrameAgainstIndex, computeFamilyCandidates, resolveParams, searchContigForMarkers,
};
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.markerGenes = exportsObj;
}
})();
