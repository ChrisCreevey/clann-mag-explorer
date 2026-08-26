const { test, report, assert } = require('./harness');
const { defaultFilters, buildBinIndex, listBinFilterOptions, applyFilters } = require('../src/model/filters');

const records = [
  { id: 'c1', length: 1000, gcContent: 0.4, codingDensity: 0.8 },
  { id: 'c2', length: 5000, gcContent: 0.6, codingDensity: 0.5 },
  { id: 'c3', length: 200, gcContent: 0.5, codingDensity: 0.9 },
];

test('defaultFilters has no active bounds and matches everything', () => {
  const result = applyFilters(records, defaultFilters());
  assert.strictEqual(result.length, 3);
});

test('length range filter excludes contigs outside the bounds', () => {
  const result = applyFilters(records, { ...defaultFilters(), lengthMin: 500, lengthMax: 4000 });
  assert.deepStrictEqual(result.map((r) => r.id), ['c1']);
});

test('GC% filter reads gcContent as a fraction but filters in percent', () => {
  const result = applyFilters(records, { ...defaultFilters(), gcMin: 55 });
  assert.deepStrictEqual(result.map((r) => r.id), ['c2']);
});

test('coding density filter works the same way', () => {
  const result = applyFilters(records, { ...defaultFilters(), codingDensityMax: 60 });
  assert.deepStrictEqual(result.map((r) => r.id), ['c2']);
});

test('buildBinIndex + bin filter selects only contigs a specific tool put in a specific bin', () => {
  const binTablesByTool = new Map([
    ['metabat2', [{ contigId: 'c1', binId: 'bin.1' }, { contigId: 'c2', binId: 'bin.2' }]],
    ['maxbin2', [{ contigId: 'c1', binId: 'binA' }]],
  ]);
  const binIndex = buildBinIndex(binTablesByTool);
  const options = listBinFilterOptions(binIndex);
  const metabatBin1 = options.find((o) => o.tool === 'metabat2' && o.binId === 'bin.1');
  const result = applyFilters(records, { ...defaultFilters(), binFilter: metabatBin1.value }, { binIndex });
  assert.deepStrictEqual(result.map((r) => r.id), ['c1']);
});

test('__unbinned__ matches contigs absent from every loaded bin table', () => {
  const binTablesByTool = new Map([['metabat2', [{ contigId: 'c1', binId: 'bin.1' }]]]);
  const binIndex = buildBinIndex(binTablesByTool);
  const result = applyFilters(records, { ...defaultFilters(), binFilter: '__unbinned__' }, { binIndex });
  assert.deepStrictEqual(result.map((r) => r.id).sort(), ['c2', 'c3']);
});

test('maxAgreementPercent excludes contigs with no recorded agreement value at all', () => {
  const agreementByContigId = new Map([['c1', 1.0], ['c2', 0.5]]);
  const result = applyFilters(records, { ...defaultFilters(), maxAgreementPercent: 80 }, { agreementByContigId });
  assert.deepStrictEqual(result.map((r) => r.id), ['c2']);
});

test('search matches contig ID or any loaded tool\'s bin ID, case-insensitively', () => {
  const binTablesByTool = new Map([['metabat2', [{ contigId: 'c1', binId: 'Foo_bin' }]]]);
  const binIndex = buildBinIndex(binTablesByTool);
  const bySearch = applyFilters(records, { ...defaultFilters(), searchText: 'foo' }, { binIndex });
  assert.deepStrictEqual(bySearch.map((r) => r.id), ['c1']);
  const byId = applyFilters(records, { ...defaultFilters(), searchText: 'C2' }, { binIndex });
  assert.deepStrictEqual(byId.map((r) => r.id), ['c2']);
});

test('filters combine with AND semantics', () => {
  const result = applyFilters(records, { ...defaultFilters(), lengthMin: 100, gcMax: 55 });
  assert.deepStrictEqual(result.map((r) => r.id).sort(), ['c1', 'c3']);
});

report();
