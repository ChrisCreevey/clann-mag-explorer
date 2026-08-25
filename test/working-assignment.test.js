const { test, report, assert } = require('./harness');
const {
  deriveInitialAssignment, reassignContigs, contigsInBin, listBinIds,
  generateNewBinId, mergeBins, assignmentToRows,
} = require('../src/model/working-assignment');

test('deriveInitialAssignment uses reconciliation majority vote when present, skipping unassigned contigs', () => {
  const reconciliation = {
    contigAgreement: [
      { contigId: 'c1', majorityMagId: 'MAG_1' },
      { contigId: 'c2', majorityMagId: 'MAG_1' },
      { contigId: 'c3', majorityMagId: null }, // no majority -> stays unassigned
    ],
  };
  const assignment = deriveInitialAssignment(new Map([['toolA', []], ['toolB', []]]), reconciliation);
  assert.strictEqual(assignment.get('c1'), 'MAG_1');
  assert.strictEqual(assignment.get('c2'), 'MAG_1');
  assert.strictEqual(assignment.has('c3'), false);
});

test('deriveInitialAssignment falls back to the single table when there is no reconciliation', () => {
  const table = [{ contigId: 'c1', binId: 'bin.1' }, { contigId: 'c2', binId: 'bin.2' }];
  const assignment = deriveInitialAssignment(new Map([['toolA', table]]), null);
  assert.strictEqual(assignment.get('c1'), 'bin.1');
  assert.strictEqual(assignment.get('c2'), 'bin.2');
});

test('deriveInitialAssignment returns an empty map when nothing was loaded', () => {
  const assignment = deriveInitialAssignment(null, null);
  assert.strictEqual(assignment.size, 0);
});

test('reassignContigs moves only the given contigs, leaves everything else untouched', () => {
  const original = new Map([['c1', 'bin.1'], ['c2', 'bin.1'], ['c3', 'bin.2']]);
  const next = reassignContigs(original, ['c1'], 'bin.2');
  assert.strictEqual(next.get('c1'), 'bin.2');
  assert.strictEqual(next.get('c2'), 'bin.1'); // untouched
  assert.strictEqual(original.get('c1'), 'bin.1'); // original map is not mutated
});

test('contigsInBin and listBinIds reflect the current assignment', () => {
  const assignment = new Map([['c1', 'bin.1'], ['c2', 'bin.1'], ['c3', 'bin.2']]);
  assert.deepStrictEqual(contigsInBin(assignment, 'bin.1').sort(), ['c1', 'c2']);
  assert.deepStrictEqual(listBinIds(assignment), ['bin.1', 'bin.2']);
});

test('generateNewBinId avoids collisions, including with previously auto-generated names', () => {
  const assignment = new Map([['c1', 'bin.new']]);
  const first = generateNewBinId(assignment);
  assert.strictEqual(first, 'bin.new.2');
  const withSecond = reassignContigs(assignment, ['c2'], first);
  const second = generateNewBinId(withSecond);
  assert.strictEqual(second, 'bin.new.3');
});

test('generateNewBinId returns the bare prefix when it is not yet used', () => {
  assert.strictEqual(generateNewBinId(new Map()), 'bin.new');
});

test('mergeBins moves every contig from the source bin into the target, target contigs untouched', () => {
  const assignment = new Map([['c1', 'bin.1'], ['c2', 'bin.1'], ['c3', 'bin.2']]);
  const merged = mergeBins(assignment, 'bin.1', 'bin.2');
  assert.strictEqual(merged.get('c1'), 'bin.2');
  assert.strictEqual(merged.get('c2'), 'bin.2');
  assert.strictEqual(merged.get('c3'), 'bin.2');
  assert.deepStrictEqual(listBinIds(merged), ['bin.2']);
});

test('assignmentToRows round-trips back to the {contigId,binId}[] shape', () => {
  const assignment = new Map([['c1', 'bin.1'], ['c2', 'bin.2']]);
  const rows = assignmentToRows(assignment);
  assert.strictEqual(rows.length, 2);
  assert.ok(rows.some((r) => r.contigId === 'c1' && r.binId === 'bin.1'));
});

report();
