(function () {
  'use strict';

// Redundancy check across finished MAGs (brief §Comparison and QC across
// the full set — "flag likely duplicate genomes split across separate
// bins"): pairwise similarity between putative MAGs (Phase 5's
// bin-reconciliation.js output), ranked so the most-likely-duplicate
// pairs surface first.
//
// Scope decision: similarity is composition-only (mean tetranucleotide
// vector + mean GC), not marker-family overlap. Unlike bin-level
// completeness/redundancy (bin-summary.js), where family *presence*
// versus *duplication* is the right question, here the families
// themselves are useless for telling organisms apart — the 40-family
// reference set is universal single-copy genes, so two MAGs from
// completely different species will both hit most of the same families
// at reasonable completeness. Composition (tetranucleotide signature)
// and GC, the same signal outliers.js already uses to separate one
// organism's contigs from another's *within* a bin, is what actually
// discriminates organism identity here.

function meanVector(vectors) {
  if (vectors.length === 0) return [];
  const dim = vectors[0].length;
  const sum = new Array(dim).fill(0);
  for (const v of vectors) for (let i = 0; i < dim; i++) sum[i] += v[i];
  return sum.map((s) => s / vectors.length);
}

function cosineSimilarity(a, b) {
  let dot = 0, normA = 0, normB = 0;
  for (let i = 0; i < a.length; i++) { dot += a[i] * b[i]; normA += a[i] * a[i]; normB += b[i] * b[i]; }
  if (normA === 0 || normB === 0) return 0;
  return dot / (Math.sqrt(normA) * Math.sqrt(normB));
}

/**
 * @param {object[]} contigRecords - per-contig records, each needs `id`,
 *   `gcContent`, and `composition` (object keyed by canonical k-mer)
 * @param {{magId:string, coreContigIds:string[], disputedContigIds:string[]}[]} putativeMags
 *   Phase 5's reconcileBins output; only MAGs with at least one contig
 *   contribute a comparison point
 * @param {{minContigs?:number, similarityThreshold?:number}} [options]
 *   minContigs: MAGs with fewer contigs than this are skipped (too little
 *   composition signal to compare meaningfully). similarityThreshold:
 *   cosine similarity above which a pair is flagged as likely-duplicate.
 * @returns {Array<{magIdA:string, magIdB:string, compositionSimilarity:number,
 *   meanGcDiff:number, likelyDuplicate:boolean}>} ranked most-similar first
 */
function computeMagRedundancy(contigRecords, putativeMags, options = {}) {
  const minContigs = options.minContigs ?? 2;
  const similarityThreshold = options.similarityThreshold ?? 0.95;

  const recordsById = new Map(contigRecords.map((r) => [r.id, r]));
  const magProfiles = [];
  for (const mag of putativeMags) {
    const contigIds = [...mag.coreContigIds, ...mag.disputedContigIds];
    const contigs = contigIds.map((id) => recordsById.get(id)).filter(Boolean);
    if (contigs.length < minContigs) continue;

    const compositionKeys = Object.keys(contigs[0].composition);
    const vectors = contigs.map((c) => compositionKeys.map((k) => c.composition[k] || 0));
    const meanGc = contigs.reduce((sum, c) => sum + c.gcContent, 0) / contigs.length;
    magProfiles.push({ magId: mag.magId, compositionVector: meanVector(vectors), meanGc });
  }

  const pairs = [];
  for (let i = 0; i < magProfiles.length; i++) {
    for (let j = i + 1; j < magProfiles.length; j++) {
      const a = magProfiles[i], b = magProfiles[j];
      const compositionSimilarity = cosineSimilarity(a.compositionVector, b.compositionVector);
      const meanGcDiff = Math.abs(a.meanGc - b.meanGc);
      pairs.push({
        magIdA: a.magId, magIdB: b.magId, compositionSimilarity, meanGcDiff,
        likelyDuplicate: compositionSimilarity >= similarityThreshold,
      });
    }
  }

  pairs.sort((x, y) => y.compositionSimilarity - x.compositionSimilarity);
  return pairs;
}

const exportsObj = { computeMagRedundancy };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.magRedundancy = exportsObj;
}
})();
