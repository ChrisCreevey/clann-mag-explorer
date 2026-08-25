(function () {
  'use strict';

// Statistical outlier detection within a bin (brief §Outlier and
// disagreement flagging): a contig sitting far from its bin's
// composition/coverage centroid. Computed per bin, independent of which
// bin table produced the grouping (works the same for a single tool's
// bins or a Phase 5 putative MAG's contig set) — callers just pass
// whichever contig list they want the centroid computed over.

function centroid(vectors) {
  if (vectors.length === 0) return [];
  const dim = vectors[0].length;
  const sum = new Array(dim).fill(0);
  for (const v of vectors) for (let i = 0; i < dim; i++) sum[i] += v[i];
  return sum.map((s) => s / vectors.length);
}

function euclideanDistance(a, b) {
  let sumSq = 0;
  for (let i = 0; i < a.length; i++) { const d = a[i] - b[i]; sumSq += d * d; }
  return Math.sqrt(sumSq);
}

/**
 * Distances from centroid, expressed as a z-score (standard deviations
 * above the bin's own mean distance) so bins with naturally tighter or
 * looser internal spread are still comparable on the same scale — a
 * contig 2 units from a tightly-clustered bin's centroid is a much
 * stronger outlier signal than the same raw distance in a naturally
 * diffuse bin.
 * @returns {Map<string, number>} contigId -> z-score (0 if the bin has
 *   only 1 contig or zero spread, since "outlier relative to itself" is undefined)
 */
function computeCentroidZScores(contigIds, vectorById) {
  const vectors = contigIds.map((id) => vectorById.get(id)).filter(Boolean);
  const c = centroid(vectors);
  const distances = contigIds.map((id) => {
    const v = vectorById.get(id);
    return v ? euclideanDistance(v, c) : null;
  });

  const validDistances = distances.filter((d) => d !== null);
  const meanDist = validDistances.reduce((s, d) => s + d, 0) / (validDistances.length || 1);
  const variance = validDistances.reduce((s, d) => s + (d - meanDist) ** 2, 0) / (validDistances.length || 1);
  const stdDev = Math.sqrt(variance);

  const zById = new Map();
  contigIds.forEach((id, i) => {
    const d = distances[i];
    zById.set(id, d === null || stdDev === 0 ? 0 : (d - meanDist) / stdDev);
  });
  return zById;
}

/**
 * @param {object[]} binContigs - per-contig records; each needs
 *   `composition` (object keyed by canonical k-mer) always, and
 *   `coverageDepths` (number[], one per sample) if a coverage table was loaded
 * @returns {Map<string, {compositionZ:number, coverageZ:number|null, combinedZ:number}>}
 */
function computeBinOutliers(binContigs) {
  const contigIds = binContigs.map((c) => c.id);

  const compositionKeys = binContigs.length ? Object.keys(binContigs[0].composition) : [];
  const compositionVectorById = new Map(
    binContigs.map((c) => [c.id, compositionKeys.map((k) => c.composition[k] || 0)])
  );
  const compositionZById = computeCentroidZScores(contigIds, compositionVectorById);

  const hasCoverage = binContigs.every((c) => Array.isArray(c.coverageDepths) && c.coverageDepths.length > 0);
  let coverageZById = null;
  if (hasCoverage) {
    const coverageVectorById = new Map(binContigs.map((c) => [c.id, c.coverageDepths]));
    coverageZById = computeCentroidZScores(contigIds, coverageVectorById);
  }

  const result = new Map();
  for (const id of contigIds) {
    const compositionZ = compositionZById.get(id) || 0;
    const coverageZ = coverageZById ? coverageZById.get(id) || 0 : null;
    // Either signal alone is enough to flag a contig as worth a look —
    // "OR", not "average", so a contig that's only unusual on one axis
    // still surfaces rather than getting diluted by the other.
    const combinedZ = coverageZ === null ? compositionZ : Math.max(compositionZ, coverageZ);
    result.set(id, { compositionZ, coverageZ, combinedZ });
  }
  return result;
}

const exportsObj = { centroid, euclideanDistance, computeCentroidZScores, computeBinOutliers };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.outliers = exportsObj;
}
})();
