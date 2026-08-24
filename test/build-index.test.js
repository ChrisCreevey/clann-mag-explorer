const { test, report, assert } = require('./harness');
const { parseHeader } = require('../build/02-index');

// Regression test: parseHeader initially read parts[1] as the taxID
// (it's the literal string "fa" — see brief §Provisioning:
// COGID.fa.<taxID>.<locus_tag>), which silently produced taxId=0 for
// every reference sequence. Caught by testing against the real generated
// data/scg40-*.bin assets, not by this unit test in isolation — added
// here so it can't regress silently again.
test('parseHeader reads the taxID from the correct segment (COGID.fa.<taxID>.<locus_tag>)', () => {
  const { family, taxId } = parseHeader('COG0012.fa.1000565.METUNv1_03812');
  assert.strictEqual(family, 'COG0012');
  assert.strictEqual(taxId, 1000565);
});

report();
