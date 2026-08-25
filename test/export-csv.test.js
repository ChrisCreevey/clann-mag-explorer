const { test, report, assert } = require('./harness');
const { assignmentToCsv, binSummaryToCsv } = require('../src/model/export-csv');

test('assignmentToCsv writes a header row plus one row per contig', () => {
  const csv = assignmentToCsv([{ contigId: 'c1', binId: 'bin.1' }, { contigId: 'c2', binId: 'bin.2' }]);
  assert.strictEqual(csv, 'contig_id,bin_id\nc1,bin.1\nc2,bin.2\n');
});

test('assignmentToCsv quotes fields containing commas', () => {
  const csv = assignmentToCsv([{ contigId: 'c1', binId: 'bin,with,commas' }]);
  assert.strictEqual(csv, 'contig_id,bin_id\nc1,"bin,with,commas"\n');
});

test('binSummaryToCsv formats percentages and leaves agreement blank when not supplied', () => {
  const summaries = [
    { binId: 'MAG_1', contigCount: 3, totalLength: 9000, n50: 4000, meanGc: 0.512, completeness: 87.5, redundancy: 2.5, mimagTier: 'high' },
  ];
  const csv = binSummaryToCsv(summaries);
  assert.strictEqual(
    csv,
    'bin_id,contig_count,total_length_bp,n50_bp,mean_gc_pct,completeness_pct,redundancy_pct,tool_agreement_fraction,mimag_tier\n' +
    'MAG_1,3,9000,4000,51.20,87.50,2.50,,high\n'
  );
});

test('binSummaryToCsv includes tool-agreement fraction when supplied', () => {
  const summaries = [
    { binId: 'MAG_1', contigCount: 3, totalLength: 9000, n50: 4000, meanGc: 0.5, completeness: 90, redundancy: 0, mimagTier: 'high' },
  ];
  const csv = binSummaryToCsv(summaries, new Map([['MAG_1', 0.8333]]));
  assert.ok(csv.includes('0.833'));
});

report();
