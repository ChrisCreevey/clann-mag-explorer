const { test, report, assert } = require('./harness');
const {
  computeN50L50, computeCompletenessRedundancy, computeMarkerContributions, computeKrakenDisagreement,
  mimagTier, computeBinSummaries,
} = require('../src/model/bin-summary');

test('computeN50L50 on a simple descending length list', () => {
  // total = 100; cumulative reaches 50 at the second contig (60 >= 50)
  const { n50, l50 } = computeN50L50([60, 30, 10], 100);
  assert.strictEqual(n50, 60);
  assert.strictEqual(l50, 1);
});

test('computeN50L50 needs more contigs when lengths are more even', () => {
  const { n50, l50 } = computeN50L50([25, 25, 25, 25], 100);
  assert.strictEqual(n50, 25);
  assert.strictEqual(l50, 2); // cumulative 25 -> 50 at the 2nd contig
});

// recallRate: 1 isolates the raw family-counting logic from the recall
// correction below (docs/scg-blast-verification.md) — with no gap to
// correct for, completeness/redundancy reduce to the plain fractions.
test('computeCompletenessRedundancy: one hit per family, one contig each -> full completeness, no redundancy', () => {
  const binContigs = [
    { id: 'c1', markerHits: [{ family: 'COG0012' }] },
    { id: 'c2', markerHits: [{ family: 'COG0016' }] },
  ];
  const { completeness, redundancy, familiesFound, rawCompleteness } = computeCompletenessRedundancy(binContigs, 1);
  assert.strictEqual(familiesFound, 2);
  assert.strictEqual(completeness, (2 / 40) * 100);
  assert.strictEqual(rawCompleteness, (2 / 40) * 100); // no gap to correct for at recallRate=1
  assert.strictEqual(redundancy, 0);
});

test('computeCompletenessRedundancy: same family on two contigs counts as redundant', () => {
  const binContigs = [
    { id: 'c1', markerHits: [{ family: 'COG0012' }] },
    { id: 'c2', markerHits: [{ family: 'COG0012' }] },
    { id: 'c3', markerHits: [{ family: 'COG0016' }] },
  ];
  const { redundancy, familiesFound } = computeCompletenessRedundancy(binContigs, 1);
  assert.strictEqual(familiesFound, 2); // COG0012, COG0016
  assert.strictEqual(redundancy, (1 / 40) * 100); // 1 of 40 expected-detectable families is duplicated
});

test('computeCompletenessRedundancy handles contigs with no marker hits at all', () => {
  const binContigs = [{ id: 'c1' }, { id: 'c2', markerHits: [] }];
  const { completeness, redundancy, familiesFound } = computeCompletenessRedundancy(binContigs, 1);
  assert.strictEqual(familiesFound, 0);
  assert.strictEqual(completeness, 0);
  assert.strictEqual(redundancy, 0); // no division by zero
});

// ---- Recall-adjusted completeness/redundancy (docs/scg-blast-verification.md) ----

test('completeness is divided by the recall rate, not read raw off the family-hit fraction', () => {
  // 20 of 40 families found is 50% raw, but at a measured 74% recall rate
  // that's evidence of a genuinely much-more-complete genome.
  const binContigs = Array.from({ length: 20 }, (_, i) => ({ id: `c${i}`, markerHits: [{ family: `FAM${i}` }] }));
  const { completeness, rawCompleteness } = computeCompletenessRedundancy(binContigs, 0.5);
  assert.strictEqual(rawCompleteness, 50);
  assert.strictEqual(completeness, 100); // 20 / (40*0.5) = 100%, capped
});

test('completeness is capped at 100% rather than exceeding it when families found beats the expected-detectable count', () => {
  const binContigs = Array.from({ length: 30 }, (_, i) => ({ id: `c${i}`, markerHits: [{ family: `FAM${i}` }] }));
  const { completeness } = computeCompletenessRedundancy(binContigs, 0.5); // expected-detectable = 20, found 30
  assert.strictEqual(completeness, 100);
});

test('redundancy is scaled by the same expected-detectable-family denominator as completeness, not by families actually found', () => {
  // Only 2 families found (a poorly-recovered bin), both duplicated —
  // "1 of 1 found families found" would swing redundancy wildly (100%)
  // off a tiny sample; scaling by the shared expected-detectable count
  // keeps it on the same, less noise-prone scale as completeness.
  const binContigs = [
    { id: 'c1', markerHits: [{ family: 'COG0012' }] },
    { id: 'c2', markerHits: [{ family: 'COG0012' }] },
  ];
  const { redundancy } = computeCompletenessRedundancy(binContigs, 0.5); // expected-detectable = 20
  assert.strictEqual(redundancy, (1 / 20) * 100);
});

test('computeCompletenessRedundancy defaults to DEFAULT_ESTIMATED_RECALL when no recallRate is given', () => {
  const { DEFAULT_ESTIMATED_RECALL } = require('../src/model/bin-summary');
  const binContigs = [{ id: 'c1', markerHits: [{ family: 'COG0012' }] }];
  const withDefault = computeCompletenessRedundancy(binContigs);
  const withExplicit = computeCompletenessRedundancy(binContigs, DEFAULT_ESTIMATED_RECALL);
  assert.strictEqual(withDefault.completeness, withExplicit.completeness);
});

test('computeMarkerContributions: a family found on only one contig is that contig\'s unique contribution', () => {
  const binContigs = [
    { id: 'c1', markerHits: [{ family: 'COG0012' }] },
    { id: 'c2', markerHits: [{ family: 'COG0016' }] },
  ];
  const contributions = computeMarkerContributions(binContigs);
  assert.deepStrictEqual(contributions.get('c1'), { uniqueFamilies: ['COG0012'], redundantFamilies: [] });
  assert.deepStrictEqual(contributions.get('c2'), { uniqueFamilies: ['COG0016'], redundantFamilies: [] });
});

test('computeMarkerContributions: a family found on multiple contigs is redundant on all of them', () => {
  const binContigs = [
    { id: 'c1', markerHits: [{ family: 'COG0012' }] },
    { id: 'c2', markerHits: [{ family: 'COG0012' }] },
  ];
  const contributions = computeMarkerContributions(binContigs);
  assert.deepStrictEqual(contributions.get('c1'), { uniqueFamilies: [], redundantFamilies: ['COG0012'] });
  assert.deepStrictEqual(contributions.get('c2'), { uniqueFamilies: [], redundantFamilies: ['COG0012'] });
});

test('computeMarkerContributions: a contig can have both unique and redundant families at once', () => {
  const binContigs = [
    { id: 'c1', markerHits: [{ family: 'COG0012' }, { family: 'COG0016' }] },
    { id: 'c2', markerHits: [{ family: 'COG0012' }] }, // shares COG0012 with c1 -> redundant on both
  ];
  const contributions = computeMarkerContributions(binContigs);
  assert.deepStrictEqual(contributions.get('c1'), { uniqueFamilies: ['COG0016'], redundantFamilies: ['COG0012'] });
});

test('computeMarkerContributions: a contig with no marker hits contributes nothing either way', () => {
  const binContigs = [{ id: 'c1', markerHits: [] }, { id: 'c2' }];
  const contributions = computeMarkerContributions(binContigs);
  assert.deepStrictEqual(contributions.get('c1'), { uniqueFamilies: [], redundantFamilies: [] });
  assert.deepStrictEqual(contributions.get('c2'), { uniqueFamilies: [], redundantFamilies: [] });
});

test('computeKrakenDisagreement flags the minority call(s), not the majority', () => {
  const binContigs = [
    { id: 'c1', krakenTaxId: 562 }, { id: 'c2', krakenTaxId: 562 }, { id: 'c3', krakenTaxId: 999 },
  ];
  const disagreement = computeKrakenDisagreement(binContigs);
  assert.strictEqual(disagreement.get('c1'), false);
  assert.strictEqual(disagreement.get('c2'), false);
  assert.strictEqual(disagreement.get('c3'), true);
});

test('computeKrakenDisagreement returns no flags when fewer than 2 contigs have a call', () => {
  const disagreement = computeKrakenDisagreement([{ id: 'c1', krakenTaxId: 562 }, { id: 'c2' }]);
  assert.strictEqual(disagreement.size, 0);
});

test('computeKrakenDisagreement ignores contigs with no call at all', () => {
  const binContigs = [{ id: 'c1', krakenTaxId: 562 }, { id: 'c2', krakenTaxId: 562 }, { id: 'c3' }];
  const disagreement = computeKrakenDisagreement(binContigs);
  assert.strictEqual(disagreement.has('c3'), false);
  assert.strictEqual(disagreement.get('c1'), false);
});

test('mimagTier: high quality needs both completeness and contamination clear', () => {
  assert.strictEqual(mimagTier(95, 2), 'high');
  assert.strictEqual(mimagTier(95, 8), 'medium'); // completeness ok, contamination too high for "high"
  assert.strictEqual(mimagTier(60, 2), 'medium');
  assert.strictEqual(mimagTier(30, 2), 'low');
  assert.strictEqual(mimagTier(60, 15), 'low'); // contamination too high even at medium completeness
});

test('mimagTier respects custom thresholds', () => {
  assert.strictEqual(mimagTier(80, 3, { highMinCompleteness: 70 }), 'high');
});

test('computeBinSummaries joins contig records to bin assignments and aggregates per bin', () => {
  const contigRecords = [
    { id: 'c1', length: 6000, gcContent: 0.5, markerHits: [{ family: 'COG0012' }] },
    { id: 'c2', length: 4000, gcContent: 0.6, markerHits: [{ family: 'COG0016' }] },
    { id: 'c3', length: 3000, gcContent: 0.4, markerHits: [] },
  ];
  const assignments = [
    { contigId: 'c1', binId: 'bin.1' },
    { contigId: 'c2', binId: 'bin.1' },
    { contigId: 'c3', binId: 'bin.2' },
  ];
  const { summaries, unmatchedContigIds } = computeBinSummaries(contigRecords, assignments);
  assert.strictEqual(unmatchedContigIds.length, 0);
  assert.strictEqual(summaries.length, 2);

  const bin1 = summaries.find((b) => b.binId === 'bin.1');
  assert.strictEqual(bin1.contigCount, 2);
  assert.strictEqual(bin1.totalLength, 10000);
  assert.ok(Math.abs(bin1.meanGc - 0.55) < 1e-9);
  assert.strictEqual(bin1.familiesFound, 2);

  // sorted by totalLength descending
  assert.strictEqual(summaries[0].binId, 'bin.1');
});

test('computeBinSummaries reports contig IDs from the bin table that have no matching loaded contig', () => {
  const contigRecords = [{ id: 'c1', length: 1000, gcContent: 0.5, markerHits: [] }];
  const assignments = [
    { contigId: 'c1', binId: 'bin.1' },
    { contigId: 'ghost_contig', binId: 'bin.1' },
  ];
  const { unmatchedContigIds } = computeBinSummaries(contigRecords, assignments);
  assert.deepStrictEqual(unmatchedContigIds, ['ghost_contig']);
});

report();
