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
