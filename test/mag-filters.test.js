const { test, report, assert } = require('./harness');
const { defaultMagFilters, applyMagFilters } = require('../src/model/mag-filters');

const mags = [
  { magId: 'MAG_1', coreCount: 10, disputedCount: 2, completeness: 92, redundancy: 1, tier: 'high', tools: ['metabat2', 'maxbin2'] },
  { magId: 'MAG_2', coreCount: 3, disputedCount: 5, completeness: 55, redundancy: 8, tier: 'medium', tools: ['metabat2'] },
  { magId: 'MAG_3', coreCount: 1, disputedCount: 0, completeness: 10, redundancy: 20, tier: 'low', tools: ['maxbin2'] },
];

test('defaultMagFilters matches every MAG', () => {
  assert.strictEqual(applyMagFilters(mags, defaultMagFilters()).length, 3);
});

test('magIdSearch matches by substring, case-insensitively', () => {
  const result = applyMagFilters(mags, { ...defaultMagFilters(), magIdSearch: 'mag_2' });
  assert.deepStrictEqual(result.map((m) => m.magId), ['MAG_2']);
});

test('tiers filters out unchecked tiers', () => {
  const filters = { ...defaultMagFilters(), tiers: { high: true, medium: false, low: false } };
  const result = applyMagFilters(mags, filters);
  assert.deepStrictEqual(result.map((m) => m.magId), ['MAG_1']);
});

test('completeness range filters by the reconciliation table\'s own completeness column', () => {
  const result = applyMagFilters(mags, { ...defaultMagFilters(), completenessMin: 50 });
  assert.deepStrictEqual(result.map((m) => m.magId).sort(), ['MAG_1', 'MAG_2']);
});

test('core/disputed count ranges work the same way', () => {
  const result = applyMagFilters(mags, { ...defaultMagFilters(), coreMin: 2, disputedMax: 3 });
  assert.deepStrictEqual(result.map((m) => m.magId), ['MAG_1']);
});

test('supportedByTool keeps only MAGs that tool contributed a bin to', () => {
  const result = applyMagFilters(mags, { ...defaultMagFilters(), supportedByTool: 'maxbin2' });
  assert.deepStrictEqual(result.map((m) => m.magId).sort(), ['MAG_1', 'MAG_3']);
});

test('filters combine with AND semantics', () => {
  const result = applyMagFilters(mags, { ...defaultMagFilters(), completenessMin: 50, supportedByTool: 'metabat2' });
  assert.deepStrictEqual(result.map((m) => m.magId).sort(), ['MAG_1', 'MAG_2']);
});

report();
