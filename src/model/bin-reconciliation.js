(function () {
  'use strict';

// Cross-tool bin matching and reconciliation (brief §Cross-tool
// reconciliation): matches equivalent bins across loaded contig->bin
// tables by contig overlap, computes per-contig agreement fraction, and
// builds the high-confidence core / disputed contig sets — the tool's
// core feature, working at contig granularity rather than DAS_Tool's
// whole-bin granularity (brief §Background).

// Contigs a tool explicitly labelled this way are treated as "this tool
// has no opinion" rather than as a real bin to match against other
// tools' bins — an "unbinned" bucket from two different tools usually
// has huge contig overlap by chance (it's just "everything left over"),
// which would otherwise make bin-matching merge them into one bogus
// putative MAG. Deliberately conservative (exact names only, not e.g.
// any purely-numeric ID) so a tool's real bin "0" isn't swept in by
// mistake — a known limitation, not exhaustive across every tool's
// convention.
const UNBINNED_LABELS = new Set(['unbinned', 'none', 'na', 'n/a', 'noclass', 'no_bin']);

function isUnbinnedLabel(binId) {
  return UNBINNED_LABELS.has(binId.trim().toLowerCase());
}

function jaccard(setA, setB) {
  let intersection = 0;
  for (const x of setA) if (setB.has(x)) intersection++;
  const union = setA.size + setB.size - intersection;
  return union === 0 ? 0 : intersection / union;
}

// Backtick+pipe separator: tool names can contain spaces, dots or
// colons (often derived from filenames), so avoid all of those.
const KEY_SEP = '|';
function binKey(tool, binId) { return `${tool}${KEY_SEP}${binId}`; }
function toolOfBinKey(k) { return k.slice(0, k.indexOf(KEY_SEP)); }

// Minimal union-find over string keys — bin counts are small (tens to a
// few hundred across all loaded tools), so no need for anything fancier.
class DisjointSet {
  constructor() { this.parent = new Map(); }
  find(x) {
    if (!this.parent.has(x)) this.parent.set(x, x);
    let root = x;
    while (this.parent.get(root) !== root) root = this.parent.get(root);
    let cur = x;
    while (this.parent.get(cur) !== root) {
      const next = this.parent.get(cur);
      this.parent.set(cur, root);
      cur = next;
    }
    return root;
  }
  union(x, y) {
    const rx = this.find(x), ry = this.find(y);
    if (rx !== ry) this.parent.set(rx, ry);
  }
}

/** @param {Map<string, {contigId:string, binId:string}[]>} assignmentsByTool */
function buildBinContigSets(assignmentsByTool) {
  const binContigs = new Map(); // binKey -> Set<contigId>
  const binMeta = new Map(); // binKey -> {tool, binId}
  for (const [tool, assignments] of assignmentsByTool) {
    for (const { contigId, binId } of assignments) {
      if (isUnbinnedLabel(binId)) continue;
      const k = binKey(tool, binId);
      if (!binContigs.has(k)) { binContigs.set(k, new Set()); binMeta.set(k, { tool, binId }); }
      binContigs.get(k).add(contigId);
    }
  }
  return { binContigs, binMeta };
}

/**
 * Mutual-best-match edges between every pair of tools: bin A (tool 1) and
 * bin B (tool 2) are linked only if each is the other's single best
 * Jaccard match AND that score clears `minJaccard` — reciprocal-best-hit
 * logic, the same conservative principle the marker-gene module uses for
 * paralog safety (brief), applied here so two only-loosely-overlapping
 * bins don't get merged into one putative MAG on a weak signal.
 */
function findMutualBestMatchEdges(binKeys, binContigs, minJaccard) {
  const toolOf = toolOfBinKey;
  const tools = [...new Set(binKeys.map(toolOf))];
  const edges = [];

  for (let i = 0; i < tools.length; i++) {
    for (let j = i + 1; j < tools.length; j++) {
      const binsA = binKeys.filter((k) => toolOf(k) === tools[i]);
      const binsB = binKeys.filter((k) => toolOf(k) === tools[j]);

      const bestForA = new Map(); // aKey -> {bKey, score}
      for (const a of binsA) {
        let best = null;
        for (const b of binsB) {
          const score = jaccard(binContigs.get(a), binContigs.get(b));
          if (score > 0 && (!best || score > best.score)) best = { bKey: b, score };
        }
        if (best) bestForA.set(a, best);
      }
      const bestForB = new Map(); // bKey -> {aKey, score}
      for (const b of binsB) {
        let best = null;
        for (const a of binsA) {
          const score = jaccard(binContigs.get(a), binContigs.get(b));
          if (score > 0 && (!best || score > best.score)) best = { aKey: a, score };
        }
        if (best) bestForB.set(b, best);
      }

      for (const [a, { bKey: b, score }] of bestForA) {
        if (score < minJaccard) continue;
        const reciprocal = bestForB.get(b);
        if (reciprocal && reciprocal.aKey === a) edges.push({ a, b, score });
      }
    }
  }
  return edges;
}

const DEFAULT_OPTIONS = {
  minJaccard: 0.1, // minimum reciprocal-best-hit overlap before two bins are considered the same putative MAG
};

/**
 * @param {Map<string, {contigId:string, binId:string}[]>} assignmentsByTool
 *   one entry per loaded binning tool, keyed by tool name
 * @param {object} [options]
 * @returns {{
 *   tools: string[],
 *   putativeMags: Array<{magId:string, members:Array<{tool:string,binId:string,contigCount:number}>,
 *     coreContigIds:string[], disputedContigIds:string[]}>,
 *   contigAgreement: Array<{contigId:string, majorityMagId:string|null, agreementFraction:number,
 *     totalVotes:number, distinctGroupsVoted:number, votes:Record<string,string|null>}>,
 *   disputedContigsRanked: object[] - contigAgreement entries with totalVotes>=2 and agreementFraction<1,
 *     most-split first
 * }}
 */
function reconcileBins(assignmentsByTool, options = {}) {
  const opts = { ...DEFAULT_OPTIONS, ...options };
  const tools = [...assignmentsByTool.keys()];
  const { binContigs, binMeta } = buildBinContigSets(assignmentsByTool);
  const binKeys = [...binContigs.keys()];

  const dsu = new DisjointSet();
  for (const k of binKeys) dsu.find(k); // register every bin, including ones with no match at all
  for (const { a, b } of findMutualBestMatchEdges(binKeys, binContigs, opts.minJaccard)) dsu.union(a, b);

  const groupsByRoot = new Map();
  for (const k of binKeys) {
    const root = dsu.find(k);
    if (!groupsByRoot.has(root)) groupsByRoot.set(root, []);
    groupsByRoot.get(root).push(k);
  }

  // Order putative MAGs by total contig count descending, so the largest
  // (most likely to be a real genome, not fragmentary agreement) leads.
  const putativeMags = [...groupsByRoot.values()]
    .map((memberKeys, i) => {
      const members = memberKeys.map((k) => ({ ...binMeta.get(k), contigCount: binContigs.get(k).size }));
      return { memberKeys, members };
    })
    .sort((a, b) => {
      const totalA = new Set(a.memberKeys.flatMap((k) => [...binContigs.get(k)])).size;
      const totalB = new Set(b.memberKeys.flatMap((k) => [...binContigs.get(k)])).size;
      return totalB - totalA;
    })
    .map((mag, i) => ({ magId: `MAG_${i + 1}`, members: mag.members, memberKeys: mag.memberKeys, coreContigIds: [], disputedContigIds: [] }));

  const magIdByBinKey = new Map();
  for (const mag of putativeMags) for (const k of mag.memberKeys) magIdByBinKey.set(k, mag.magId);

  // Per-contig votes: which putative MAG each tool assigned this contig
  // to, skipping tools that left it unbinned/didn't mention it at all
  // (both mean "no opinion", not "voted for no MAG").
  const contigVotes = new Map(); // contigId -> Map<tool, magId>
  for (const [tool, assignments] of assignmentsByTool) {
    for (const { contigId, binId } of assignments) {
      if (isUnbinnedLabel(binId)) continue;
      const magId = magIdByBinKey.get(binKey(tool, binId));
      if (!contigVotes.has(contigId)) contigVotes.set(contigId, new Map());
      contigVotes.get(contigId).set(tool, magId);
    }
  }

  const magsById = new Map(putativeMags.map((m) => [m.magId, m]));
  const contigAgreement = [];
  for (const [contigId, votes] of contigVotes) {
    const counts = new Map();
    for (const magId of votes.values()) counts.set(magId, (counts.get(magId) || 0) + 1);
    let majorityMagId = null, majorityCount = 0;
    for (const [magId, count] of counts) {
      if (count > majorityCount) { majorityCount = count; majorityMagId = magId; }
    }
    const totalVotes = votes.size;
    const agreementFraction = totalVotes ? majorityCount / totalVotes : 0;
    contigAgreement.push({
      contigId, majorityMagId, agreementFraction, totalVotes,
      distinctGroupsVoted: counts.size, votes: Object.fromEntries(votes),
    });

    if (majorityMagId) {
      const mag = magsById.get(majorityMagId);
      if (agreementFraction === 1) mag.coreContigIds.push(contigId);
      else mag.disputedContigIds.push(contigId);
    }
  }

  const disputedContigsRanked = contigAgreement
    .filter((c) => c.totalVotes >= 2 && c.agreementFraction < 1)
    .sort((a, b) => a.agreementFraction - b.agreementFraction || b.distinctGroupsVoted - a.distinctGroupsVoted);

  for (const mag of putativeMags) delete mag.memberKeys; // internal-only, not part of the public shape

  return { tools, putativeMags, contigAgreement, disputedContigsRanked };
}

const exportsObj = { reconcileBins, isUnbinnedLabel, jaccard };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.binReconciliation = exportsObj;
}
})();
