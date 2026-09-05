(function () {
  'use strict';

// Best-attempt contig ID matching between a loaded contig->bin table and
// the actual assembly: some binning tools rewrite contig IDs before
// clustering (CONCOCT's default `cut_up_fasta.py` step splits every
// contig into fixed-length chunks and clusters those, producing IDs like
// `<contig>.concoct_part_0`, `.concoct_part_1`, ... instead of the
// original contig ID), and some write their columns in the opposite
// order this app assumes (VAMB's native cluster output is
// clustername\tcontigname, not this app's assumed contig\tbin) — both of
// which otherwise make every one of that tool's rows silently fail to
// match the assembly, not because it was run against a different
// assembly, but because of a well-known, mechanical format quirk. This
// module detects both cases (and is structured to add more known
// conventions later) and, only when a candidate fix recovers a high
// match rate, rewrites the table's contig IDs and collapses any contig
// that split into multiple rows back into one (majority vote across its
// parts' bin calls, ties broken by first-seen order for determinism).
//
// If NO known convention gets the match rate high, that's the genuine
// "this table doesn't belong to this assembly" signal the brief's
// "different assembly version" concern is about — surfaced as a report
// rather than guessed at further, since inventing more pattern-matching
// heuristics for an actually-mismatched assembly would just produce
// wrong-but-confident results.

const KNOWN_SUFFIX_PATTERNS = [
  { name: 'concoct_part', regex: /\.concoct_part_\d+$/, label: 'CONCOCT contig-splitting suffix (.concoct_part_N)' },
];

const HIGH_MATCH_THRESHOLD = 0.9; // fraction of distinct contig IDs that must resolve against the assembly to trust a strategy

function matchRate(ids, referenceIds) {
  if (ids.length === 0) return 0;
  let matched = 0;
  for (const id of ids) if (referenceIds.has(id)) matched++;
  return matched / ids.length;
}

/** Collapses possibly-duplicate {contigId, binId} rows (after stripping a suffix) via majority vote per contigId. */
function collapseByMajorityVote(assignments) {
  const votesByContig = new Map(); // contigId -> Map<binId, count>
  const firstSeenOrder = new Map(); // contigId -> binId[] in first-seen order, for deterministic tie-break
  for (const { contigId, binId } of assignments) {
    if (!votesByContig.has(contigId)) { votesByContig.set(contigId, new Map()); firstSeenOrder.set(contigId, []); }
    const votes = votesByContig.get(contigId);
    if (!votes.has(binId)) firstSeenOrder.get(contigId).push(binId);
    votes.set(binId, (votes.get(binId) || 0) + 1);
  }

  const result = [];
  let collapsedCount = 0;
  for (const [contigId, votes] of votesByContig) {
    if (votes.size > 1) collapsedCount++;
    let bestBinId = null, bestCount = -1;
    for (const binId of firstSeenOrder.get(contigId)) {
      const count = votes.get(binId);
      if (count > bestCount) { bestCount = count; bestBinId = binId; }
    }
    result.push({ contigId, binId: bestBinId });
  }
  return { assignments: result, collapsedCount };
}

/**
 * @param {{contigId:string, binId:string}[]} assignments - one tool's loaded table
 * @param {Set<string>} referenceContigIds - every contig ID the loaded assembly actually has
 * @returns {{assignments: object[], report: {
 *   strategy: 'exact'|string, applied: boolean,
 *   matchRateBefore: number, matchRateAfter: number,
 *   collapsedCount: number, patternLabel: string|null,
 * }}}
 */
function bestAttemptRemapAssignments(assignments, referenceContigIds) {
  const originalIds = assignments.map((a) => a.contigId);
  const matchRateBefore = matchRate(originalIds, referenceContigIds);

  if (matchRateBefore >= HIGH_MATCH_THRESHOLD) {
    return {
      assignments,
      report: { strategy: 'exact', applied: false, matchRateBefore, matchRateAfter: matchRateBefore, collapsedCount: 0, patternLabel: null },
    };
  }

  // Each candidate below is {strategy, rate, patternLabel, buildAssignments()}
  // — same shape regardless of which recovery it represents, so they can
  // all compete on `rate` and the best one wins, same principle as picking
  // the best-scoring reciprocal match elsewhere in this app.
  let best = null;

  for (const pattern of KNOWN_SUFFIX_PATTERNS) {
    const strippedIds = originalIds.map((id) => id.replace(pattern.regex, ''));
    const rate = matchRate(strippedIds, referenceContigIds);
    if (rate > matchRateBefore && (!best || rate > best.rate)) {
      best = {
        strategy: pattern.name, rate, patternLabel: pattern.label,
        buildAssignments: () => assignments.map((a, i) => ({ contigId: strippedIds[i], binId: a.binId })),
      };
    }
  }

  // Some tools write bin\tcontig rather than this app's assumed
  // contig\tbin — VAMB's native cluster output (clustername, contigname),
  // notably, the exact reverse of DAS_Tool's Fasta_to_Contig2Bin.sh
  // convention everything here is otherwise built around. Parsed the
  // wrong way round, every real contig ID lands in the *binId* field, so
  // each becomes its own one-contig "bin" that can never reciprocal-match
  // anything — a real, confusing symptom this shipped with: that tool
  // showing zero contribution to every putative MAG, and every contig it
  // touches counted "disputed" purely because its phantom singleton bins
  // never agree with anyone. Detected the same content-based way as every
  // other strategy here: if reading the *second* column as the contig ID
  // recovers a high match rate against the real assembly, the columns are
  // almost certainly swapped.
  const swappedIds = assignments.map((a) => a.binId);
  const swappedRate = matchRate(swappedIds, referenceContigIds);
  if (swappedRate > matchRateBefore && (!best || swappedRate > best.rate)) {
    best = {
      strategy: 'swapped-columns', rate: swappedRate,
      patternLabel: 'bin→contig column order (e.g. VAMB’s native cluster output), not contig→bin',
      buildAssignments: () => assignments.map((a) => ({ contigId: a.binId, binId: a.contigId })),
    };
  }

  if (best && best.rate >= HIGH_MATCH_THRESHOLD) {
    const { assignments: collapsed, collapsedCount } = collapseByMajorityVote(best.buildAssignments());
    return {
      assignments: collapsed,
      report: {
        strategy: best.strategy, applied: true,
        matchRateBefore, matchRateAfter: best.rate, collapsedCount, patternLabel: best.patternLabel,
      },
    };
  }

  // Nothing recovered a high match rate — leave the table as loaded and
  // report the low match rate as-is, rather than apply a weak heuristic.
  return {
    assignments,
    report: {
      strategy: 'none', applied: false,
      matchRateBefore, matchRateAfter: best ? best.rate : matchRateBefore, collapsedCount: 0, patternLabel: null,
    },
  };
}

const exportsObj = { matchRate, bestAttemptRemapAssignments, KNOWN_SUFFIX_PATTERNS, HIGH_MATCH_THRESHOLD };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.contigIdMatching = exportsObj;
}
})();
