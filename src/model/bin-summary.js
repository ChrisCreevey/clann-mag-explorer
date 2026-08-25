(function () {
  'use strict';

// Per-bin summary statistics (brief §Per-bin / per-putative-MAG summary):
// contig count, total length, N50/L50, mean GC, SCG-based
// completeness/redundancy from Phase 3's marker-gene tags, and a
// MIMAG-style quality tier — the first thing built once a single
// contig->bin table is loaded (brief's Phase 4).

const TOTAL_MARKER_FAMILIES = 40; // fixed by the shipped reference set (brief §Provisioning)

// MIMAG-style tiers (Bowers et al. 2017), completeness/contamination only —
// this module has no rRNA/tRNA signal to check the real MIMAG standard's
// full criteria, so these are a simplified proxy, not a certification.
// Exposed as parameters (not hardcoded) per the brief's "thresholds shown
// and adjustable" framing, matching DEFAULT_PARAMS-as-overridable-object
// used elsewhere (e.g. marker-genes.js).
const DEFAULT_MIMAG_THRESHOLDS = {
  highMinCompleteness: 90, highMaxContamination: 5,
  mediumMinCompleteness: 50, mediumMaxContamination: 10,
};

function computeN50L50(lengthsDesc, totalLength) {
  let cumulative = 0;
  for (let i = 0; i < lengthsDesc.length; i++) {
    cumulative += lengthsDesc[i];
    if (cumulative >= totalLength / 2) return { n50: lengthsDesc[i], l50: i + 1 };
  }
  return { n50: 0, l50: 0 };
}

/**
 * Completeness = fraction of the 40 marker families found anywhere in the
 * bin. Redundancy = fraction of the families that were found (not of all
 * 40) that showed up on more than one contig — the standard CheckM-style
 * contamination proxy (extra copies beyond the expected single copy per
 * family), computed here from Phase 3's per-contig calls rather than a
 * separate whole-bin re-search, since brief §Marker-gene identification
 * module frames per-contig tags as the thing everything downstream
 * aggregates over.
 * @param {object[]} binContigs - this bin's per-contig records, each with
 *   an optional `markerHits: {family: string}[]`
 */
function computeCompletenessRedundancy(binContigs) {
  const contigsByFamily = new Map(); // family -> Set<contigId>
  for (const contig of binContigs) {
    for (const hit of contig.markerHits || []) {
      if (!contigsByFamily.has(hit.family)) contigsByFamily.set(hit.family, new Set());
      contigsByFamily.get(hit.family).add(contig.id);
    }
  }
  const familiesFound = contigsByFamily.size;
  let familiesWithMultipleContigs = 0;
  for (const contigSet of contigsByFamily.values()) {
    if (contigSet.size > 1) familiesWithMultipleContigs++;
  }

  const completeness = (familiesFound / TOTAL_MARKER_FAMILIES) * 100;
  const redundancy = familiesFound ? (familiesWithMultipleContigs / familiesFound) * 100 : 0;
  return { completeness, redundancy, familiesFound };
}

/**
 * Per-contig marker contribution (brief §Marker-gene identification
 * module — "the distinctive output beyond a bin-level number"): for each
 * contig, which of its called families are found on no other contig in
 * the same bin (unique — removing this contig loses completeness) versus
 * also found elsewhere in the bin (redundant — removing this contig
 * reduces apparent contamination without costing completeness). A contig
 * with high redundant and zero unique contribution is a low-risk removal
 * candidate.
 * @param {object[]} binContigs - this bin's per-contig records
 * @returns {Map<string, {uniqueFamilies:string[], redundantFamilies:string[]}>} keyed by contig id
 */
function computeMarkerContributions(binContigs) {
  const contigsByFamily = new Map(); // family -> Set<contigId>
  for (const contig of binContigs) {
    for (const hit of contig.markerHits || []) {
      if (!contigsByFamily.has(hit.family)) contigsByFamily.set(hit.family, new Set());
      contigsByFamily.get(hit.family).add(contig.id);
    }
  }

  const contributions = new Map();
  for (const contig of binContigs) {
    const uniqueFamilies = [];
    const redundantFamilies = [];
    for (const hit of contig.markerHits || []) {
      const contigsWithThisFamily = contigsByFamily.get(hit.family);
      if (contigsWithThisFamily.size > 1) redundantFamilies.push(hit.family);
      else uniqueFamilies.push(hit.family);
    }
    contributions.set(contig.id, { uniqueFamilies, redundantFamilies });
  }
  return contributions;
}

/**
 * Taxonomic-disagreement flag (brief §Outlier and disagreement flagging):
 * a contig whose Kraken2 call differs from the rest of its bin. Scope
 * decision: exact-taxID majority-vote mismatch, not lineage-aware (e.g.
 * two different species in the same genus would still count as a
 * "disagreement" here). A full lineage-aware version needs a taxonomy
 * tree covering whatever arbitrary taxIDs Kraken2's default database
 * calls — unlike the marker-gene provenance check (marker-taxonomy.js),
 * which only ever needs the ~5,000 taxa this app's own fixed 40-family
 * reference set touches, Kraken2 can call anything in NCBI's full
 * taxonomy. Shipping that whole tree is out of scope here; a real
 * lineage-aware upgrade is possible if a student also loads a full-run
 * .breport alongside the per-contig calls (breport.js already builds a
 * real taxonomy-tree.js from one, per the eDNA Explorer port) — noted as
 * a follow-up, not built now.
 * @param {object[]} binContigs - per-contig records with an optional
 *   `krakenTaxId: number` (from kraken2-contigs.js)
 * @returns {Map<string, boolean>} contigId -> true if this contig's call
 *   disagrees with the bin's majority call (only set for contigs that
 *   have a call at all, in a bin where at least 2 contigs do)
 */
function computeKrakenDisagreement(binContigs) {
  const withCalls = binContigs.filter((c) => c.krakenTaxId != null);
  const disagreement = new Map();
  if (withCalls.length < 2) return disagreement;

  const counts = new Map();
  for (const c of withCalls) counts.set(c.krakenTaxId, (counts.get(c.krakenTaxId) || 0) + 1);
  let majorityTaxId = null, majorityCount = 0;
  for (const [taxId, count] of counts) if (count > majorityCount) { majorityCount = count; majorityTaxId = taxId; }

  for (const c of withCalls) disagreement.set(c.id, c.krakenTaxId !== majorityTaxId);
  return disagreement;
}

function mimagTier(completeness, contamination, thresholds) {
  const t = { ...DEFAULT_MIMAG_THRESHOLDS, ...thresholds };
  if (completeness > t.highMinCompleteness && contamination < t.highMaxContamination) return 'high';
  if (completeness >= t.mediumMinCompleteness && contamination < t.mediumMaxContamination) return 'medium';
  return 'low';
}

/**
 * @param {object[]} contigRecords - every loaded contig's stats record
 *   (id, length, gcContent, markerHits, ...), keyed by `id`
 * @param {{contigId: string, binId: string}[]} assignments - one binning
 *   tool's contig->bin table
 * @param {{completenessSource?: 'builtin'|'supplied', thresholds?: object}} [options]
 *   completenessSource is currently always 'builtin' (the brief's
 *   "supplied pre-computed hits, preferred when present" path is Phase 4+
 *   scope not yet wired to an input parser — see docs/phase1-investigation.md)
 * @returns {object[]} one summary per bin, sorted by totalLength descending
 */
function computeBinSummaries(contigRecords, assignments, options = {}) {
  const recordsById = new Map(contigRecords.map((r) => [r.id, r]));
  const binIdToContigs = new Map();
  const unmatchedContigIds = [];

  for (const { contigId, binId } of assignments) {
    const record = recordsById.get(contigId);
    if (!record) { unmatchedContigIds.push(contigId); continue; }
    if (!binIdToContigs.has(binId)) binIdToContigs.set(binId, []);
    binIdToContigs.get(binId).push(record);
  }

  const summaries = [];
  for (const [binId, contigs] of binIdToContigs) {
    const lengthsDesc = contigs.map((c) => c.length).sort((a, b) => b - a);
    const totalLength = lengthsDesc.reduce((sum, len) => sum + len, 0);
    const { n50, l50 } = computeN50L50(lengthsDesc, totalLength);
    const meanGc = contigs.reduce((sum, c) => sum + c.gcContent, 0) / contigs.length;

    const { completeness, redundancy, familiesFound } = computeCompletenessRedundancy(contigs);
    const tier = mimagTier(completeness, redundancy, options.thresholds);

    summaries.push({
      binId, contigCount: contigs.length, totalLength, n50, l50, meanGc,
      completeness, redundancy, familiesFound, mimagTier: tier,
    });
  }

  summaries.sort((a, b) => b.totalLength - a.totalLength);
  return { summaries, unmatchedContigIds };
}

const exportsObj = {
  TOTAL_MARKER_FAMILIES, DEFAULT_MIMAG_THRESHOLDS,
  computeN50L50, computeCompletenessRedundancy, computeMarkerContributions, computeKrakenDisagreement,
  mimagTier, computeBinSummaries,
};
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.binSummary = exportsObj;
}
})();
