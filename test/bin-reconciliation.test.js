const { test, report, assert } = require('./harness');
const { reconcileBins, isUnbinnedLabel, jaccard } = require('../src/model/bin-reconciliation');

function assignments(pairs) {
  return pairs.map(([contigId, binId]) => ({ contigId, binId }));
}

test('jaccard of identical sets is 1, disjoint sets is 0', () => {
  assert.strictEqual(jaccard(new Set(['a', 'b']), new Set(['a', 'b'])), 1);
  assert.strictEqual(jaccard(new Set(['a']), new Set(['b'])), 0);
});

test('isUnbinnedLabel matches common conventions case-insensitively, not arbitrary numeric IDs', () => {
  assert.ok(isUnbinnedLabel('unbinned'));
  assert.ok(isUnbinnedLabel('Unbinned'));
  assert.ok(isUnbinnedLabel('NA'));
  assert.ok(!isUnbinnedLabel('0'));
  assert.ok(!isUnbinnedLabel('bin.1'));
});

test('two tools with identical bins produce one putative MAG with full core, no disputed', () => {
  const toolA = assignments([['c1', 'bin.1'], ['c2', 'bin.1'], ['c3', 'bin.2']]);
  const toolB = assignments([['c1', 'A'], ['c2', 'A'], ['c3', 'B']]);
  const result = reconcileBins(new Map([['toolA', toolA], ['toolB', toolB]]));

  assert.strictEqual(result.putativeMags.length, 2);
  const magWithC1 = result.putativeMags.find((m) => m.coreContigIds.includes('c1'));
  assert.ok(magWithC1, 'c1 should be in some MAG\'s core set');
  assert.ok(magWithC1.coreContigIds.includes('c2'));
  assert.strictEqual(magWithC1.disputedContigIds.length, 0);
  assert.strictEqual(magWithC1.members.length, 2); // one bin from each tool
});

test('a bin that reciprocal-matches nothing and never wins a contig\'s majority vote is dropped from putativeMags entirely', () => {
  // Real messy multi-tool input (a near-singleton binner like VAMB
  // especially) produces exactly this: two tools agree on a large bin,
  // a third tool puts one of that bin's contigs in its own tiny bin. The
  // tiny bin's overlap with the big one (1/15) falls below minJaccard, so
  // it never reciprocal-matches, and it loses every majority vote it
  // could ever cast (2 tools vs. 1) — it should end up in nobody's
  // core/disputed set, and NOT appear in putativeMags as a phantom
  // "putative MAG" with zero contigs, even though its bin ID still shows
  // up as a minority vote (contigAgreement's `votes`).
  const bigContigs = Array.from({ length: 15 }, (_, i) => `c${i + 1}`);
  const toolA = assignments(bigContigs.map((c) => [c, 'bigBin']));
  const toolB = assignments(bigContigs.map((c) => [c, 'bigBin']));
  const toolC = assignments([['c1', 'lone1']]); // only opines on c1, in its own tiny bin

  const result = reconcileBins(new Map([['toolA', toolA], ['toolB', toolB], ['toolC', toolC]]));

  assert.strictEqual(result.putativeMags.length, 1, 'lone1 should not survive as its own putative MAG');
  const bigMag = result.putativeMags[0];
  // c1 isn't unanimous (toolC dissents), so it's disputed rather than
  // core — but it still belongs to the 2-tool majority MAG, not to lone1.
  assert.ok(bigMag.disputedContigIds.includes('c1'), 'c1 goes to the 2-tool majority MAG (as disputed), not the lone tiny bin');

  const c1Agreement = result.contigAgreement.find((c) => c.contigId === 'c1');
  assert.ok(c1Agreement.votes.toolC, 'toolC\'s vote is still recorded even though its MAG was dropped from putativeMags');
  assert.notStrictEqual(c1Agreement.votes.toolC, bigMag.magId, 'toolC voted for a different (dropped) MAG than the majority');
});

test('a contig two tools disagree on is disputed, not core, and ranked in disputedContigsRanked', () => {
  const toolA = assignments([['c1', 'bin.1'], ['c2', 'bin.1'], ['c3', 'bin.1'], ['c4', 'bin.2'], ['c5', 'bin.2']]);
  const toolB = assignments([['c1', 'A'], ['c2', 'A'], ['c3', 'B'], ['c4', 'B'], ['c5', 'B']]);
  const result = reconcileBins(new Map([['toolA', toolA], ['toolB', toolB]]));

  const c3 = result.contigAgreement.find((c) => c.contigId === 'c3');
  assert.ok(c3, 'c3 should have an agreement record');
  assert.strictEqual(c3.totalVotes, 2);
  assert.ok(c3.agreementFraction < 1, 'c3 should be disputed (tools disagree on its group)');

  const ranked = result.disputedContigsRanked.map((c) => c.contigId);
  assert.ok(ranked.includes('c3'));
});

test('unbinned label from one tool contributes no vote for that contig', () => {
  const toolA = assignments([['c1', 'bin.1'], ['c2', 'bin.1']]);
  const toolB = assignments([['c1', 'unbinned'], ['c2', 'A']]);
  const result = reconcileBins(new Map([['toolA', toolA], ['toolB', toolB]]));

  const c1 = result.contigAgreement.find((c) => c.contigId === 'c1');
  assert.strictEqual(c1.totalVotes, 1); // only toolA voted
  assert.strictEqual(c1.agreementFraction, 1); // trivially "unanimous" among the tools that did vote
});

test('a bin only one tool supports (no reciprocal match) becomes its own singleton putative MAG', () => {
  const toolA = assignments([['c1', 'bin.1'], ['c2', 'bin.1']]);
  const toolB = assignments([['c3', 'X'], ['c4', 'X']]); // completely disjoint contig set from bin.1
  const result = reconcileBins(new Map([['toolA', toolA], ['toolB', toolB]]));

  assert.strictEqual(result.putativeMags.length, 2);
  for (const mag of result.putativeMags) assert.strictEqual(mag.members.length, 1);
});

test('three tools: majority-of-two-vs-one is still a disputed contig, not core', () => {
  // toolA and toolB agree bin.1/A = {c1, c2} exactly (jaccard 1 -> merges).
  // toolC agrees on c2 (its bin A2={c2} overlaps bin.1 at jaccard 0.5, still
  // merges into the same group), but puts c1 in a bin (X) padded with 9
  // contigs no other tool mentions — jaccard(X, bin.1) = 1/11 ≈ 0.09, just
  // under the default 0.1 threshold, so X stays a separate putative MAG
  // rather than getting transitively merged in on weak evidence.
  const toolA = assignments([['c1', 'bin.1'], ['c2', 'bin.1']]);
  const toolB = assignments([['c1', 'A'], ['c2', 'A']]);
  const junk = ['z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'z9'];
  const toolC = assignments([['c2', 'A2'], ...junk.map((j) => [j, 'X']), ['c1', 'X']]);
  const result = reconcileBins(new Map([['toolA', toolA], ['toolB', toolB], ['toolC', toolC]]));

  const c1 = result.contigAgreement.find((c) => c.contigId === 'c1');
  const c2 = result.contigAgreement.find((c) => c.contigId === 'c2');
  assert.strictEqual(c1.totalVotes, 3);
  assert.ok(c1.agreementFraction < 1, `c1: 2 of 3 tools agree, not unanimous (got ${c1.agreementFraction})`);
  assert.strictEqual(c2.agreementFraction, 1, 'c2: all 3 tools agree');
});

test('minJaccard option controls how loose an overlap still counts as a match', () => {
  // bin.1 (toolA) and A (toolB) share only 1 of 4 contigs each -> low overlap
  const toolA = assignments([['c1', 'bin.1'], ['c2', 'bin.1'], ['c3', 'bin.1'], ['c4', 'bin.1']]);
  const toolB = assignments([['c1', 'A'], ['c5', 'A'], ['c6', 'A'], ['c7', 'A']]);
  const strict = reconcileBins(new Map([['toolA', toolA], ['toolB', toolB]]), { minJaccard: 0.5 });
  const loose = reconcileBins(new Map([['toolA', toolA], ['toolB', toolB]]), { minJaccard: 0.05 });

  assert.strictEqual(strict.putativeMags.length, 2, 'weak overlap should not merge under a strict threshold');
  assert.strictEqual(loose.putativeMags.length, 1, 'the same overlap should merge under a loose threshold');
});

report();
