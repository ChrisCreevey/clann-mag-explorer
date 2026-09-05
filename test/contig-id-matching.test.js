const { test, report, assert } = require('./harness');
const { matchRate, bestAttemptRemapAssignments } = require('../src/model/contig-id-matching');

test('matchRate is the fraction of ids present in the reference set', () => {
  const ref = new Set(['a', 'b', 'c']);
  assert.strictEqual(matchRate(['a', 'b', 'x'], ref), 2 / 3);
  assert.strictEqual(matchRate([], ref), 0);
});

test('bestAttemptRemapAssignments leaves an already-matching table untouched', () => {
  const assignments = [{ contigId: 'a', binId: 'bin1' }, { contigId: 'b', binId: 'bin1' }];
  const ref = new Set(['a', 'b', 'c']);
  const { assignments: result, report } = bestAttemptRemapAssignments(assignments, ref);
  assert.strictEqual(report.strategy, 'exact');
  assert.strictEqual(report.applied, false);
  assert.deepStrictEqual(result, assignments);
});

test('bestAttemptRemapAssignments strips CONCOCT part suffixes when that recovers a high match rate', () => {
  const assignments = [
    { contigId: 'contigA.concoct_part_0', binId: '90' },
    { contigId: 'contigB.concoct_part_0', binId: '12' },
    { contigId: 'contigC.concoct_part_0', binId: '12' },
  ];
  const ref = new Set(['contigA', 'contigB', 'contigC']);
  const { assignments: result, report } = bestAttemptRemapAssignments(assignments, ref);
  assert.strictEqual(report.strategy, 'concoct_part');
  assert.strictEqual(report.applied, true);
  assert.strictEqual(report.matchRateAfter, 1);
  assert.deepStrictEqual(
    result.sort((a, b) => a.contigId.localeCompare(b.contigId)),
    [{ contigId: 'contigA', binId: '90' }, { contigId: 'contigB', binId: '12' }, { contigId: 'contigC', binId: '12' }]
  );
});

test('a contig split into multiple CONCOCT parts collapses to one row via majority vote', () => {
  const assignments = [
    { contigId: 'contigA.concoct_part_0', binId: '90' },
    { contigId: 'contigA.concoct_part_1', binId: '90' },
    { contigId: 'contigA.concoct_part_2', binId: '49' }, // dissenting part
  ];
  const ref = new Set(['contigA']);
  const { assignments: result, report } = bestAttemptRemapAssignments(assignments, ref);
  assert.strictEqual(result.length, 1);
  assert.strictEqual(result[0].contigId, 'contigA');
  assert.strictEqual(result[0].binId, '90'); // majority of the 3 parts
  assert.strictEqual(report.collapsedCount, 1);
});

test('bestAttemptRemapAssignments detects a bin<->contig column swap (e.g. VAMB\'s native clustername\\tcontigname output) and corrects it', () => {
  // Written the wrong way round for this app's assumed contig,bin order:
  // real bin names in the first column, real contig IDs in the second.
  const assignments = [
    { contigId: 'cluster_1', binId: 'contigA' },
    { contigId: 'cluster_1', binId: 'contigB' },
    { contigId: 'cluster_2', binId: 'contigC' },
  ];
  const ref = new Set(['contigA', 'contigB', 'contigC']);
  const { assignments: result, report } = bestAttemptRemapAssignments(assignments, ref);
  assert.strictEqual(report.strategy, 'swapped-columns');
  assert.strictEqual(report.applied, true);
  assert.strictEqual(report.matchRateAfter, 1);
  assert.deepStrictEqual(
    result.sort((a, b) => a.contigId.localeCompare(b.contigId)),
    [{ contigId: 'contigA', binId: 'cluster_1' }, { contigId: 'contigB', binId: 'cluster_1' }, { contigId: 'contigC', binId: 'cluster_2' }]
  );
});

test('a swap-detected table with a genuinely duplicated contig row still collapses via majority vote', () => {
  const assignments = [
    { contigId: 'cluster_1', binId: 'contigA' },
    { contigId: 'cluster_2', binId: 'contigA' }, // same real contig, second row disagrees (rare, but exercise the path)
  ];
  const ref = new Set(['contigA']);
  const { assignments: result, report } = bestAttemptRemapAssignments(assignments, ref);
  assert.strictEqual(report.strategy, 'swapped-columns');
  assert.strictEqual(result.length, 1);
  assert.strictEqual(result[0].contigId, 'contigA');
  assert.strictEqual(report.collapsedCount, 1);
});

test('bestAttemptRemapAssignments picks whichever strategy recovers the higher match rate when both could apply', () => {
  // Swapping alone gets every contig ID right; nothing here looks like a
  // CONCOCT suffix, so the swap strategy should win outright.
  const assignments = [
    { contigId: 'cluster_1', binId: 'contigA' },
    { contigId: 'cluster_1', binId: 'contigB' },
  ];
  const ref = new Set(['contigA', 'contigB']);
  const { report } = bestAttemptRemapAssignments(assignments, ref);
  assert.strictEqual(report.strategy, 'swapped-columns');
  assert.strictEqual(report.matchRateAfter, 1);
});

test('bestAttemptRemapAssignments reports strategy "none" and leaves the table alone when nothing recovers a high match rate', () => {
  const assignments = [
    { contigId: 'totallyDifferentContig1', binId: 'bin1' },
    { contigId: 'totallyDifferentContig2', binId: 'bin1' },
  ];
  const ref = new Set(['contigA', 'contigB']);
  const { assignments: result, report } = bestAttemptRemapAssignments(assignments, ref);
  assert.strictEqual(report.strategy, 'none');
  assert.strictEqual(report.applied, false);
  assert.strictEqual(report.matchRateAfter, 0);
  assert.deepStrictEqual(result, assignments);
});

report();
