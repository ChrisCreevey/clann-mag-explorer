const { test, report, assert } = require('./harness');
const { pickComparisonBins, buildMagQcPoints } = require('../src/model/qc-comparison');

test('pickComparisonBins picks highest and lowest completeness-minus-redundancy score', () => {
  const summaries = [
    { binId: 'good', contigCount: 5, completeness: 95, redundancy: 2 },
    { binId: 'mid', contigCount: 3, completeness: 60, redundancy: 20 },
    { binId: 'bad', contigCount: 4, completeness: 20, redundancy: 40 },
  ];
  const { good, bad } = pickComparisonBins(summaries);
  assert.strictEqual(good.binId, 'good');
  assert.strictEqual(bad.binId, 'bad');
});

test('pickComparisonBins excludes single-contig bins', () => {
  const summaries = [
    { binId: 'solo', contigCount: 1, completeness: 100, redundancy: 0 },
    { binId: 'ok', contigCount: 3, completeness: 50, redundancy: 10 },
  ];
  const { good, bad } = pickComparisonBins(summaries);
  assert.strictEqual(good, null);
  assert.strictEqual(bad, null);
});

test('pickComparisonBins returns null when fewer than 2 eligible bins', () => {
  const summaries = [{ binId: 'only', contigCount: 5, completeness: 90, redundancy: 1 }];
  const { good, bad } = pickComparisonBins(summaries);
  assert.strictEqual(good, null);
  assert.strictEqual(bad, null);
});

test('buildMagQcPoints maps summaries to completeness/redundancy points coloured by supporting tools', () => {
  const magSummaries = [
    { binId: 'MAG_1', completeness: 90, redundancy: 5 },
    { binId: 'MAG_2', completeness: 40, redundancy: 20 },
  ];
  const putativeMags = [
    { magId: 'MAG_1', members: [{ tool: 'metabat2', binId: '3' }, { tool: 'maxbin2', binId: '7' }] },
    { magId: 'MAG_2', members: [{ tool: 'metabat2', binId: '9' }] },
  ];
  const points = buildMagQcPoints(magSummaries, putativeMags);
  assert.deepStrictEqual(points, [
    { id: 'MAG_1', x: 90, y: 5, colorKey: 'maxbin2, metabat2' },
    { id: 'MAG_2', x: 40, y: 20, colorKey: 'metabat2' },
  ]);
});

report();
