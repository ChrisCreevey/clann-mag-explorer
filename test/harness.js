// Minimal zero-dependency test harness, matching the suite's "no external
// dependencies" convention. Usage in a test file:
//
//   const { test, report } = require('./harness');
//   test('does the thing', () => { assert.strictEqual(1 + 1, 2); });
//   test('does the async thing', async () => { assert.ok(await somePromise); });
//   report();
//
// Tests run in registration order via a single chained promise shared
// across every file `require`-ing this module (Node caches the module, so
// `pending`/`passed`/`failed` are one instance for the whole `test/run.js`
// run) — a sync test's body just resolves that link immediately, so mixing
// sync and async test() calls in the same or different files stays ordered
// and doesn't need special-casing.

const assert = require('assert');

let passed = 0;
let failed = 0;
const failures = [];
let pending = Promise.resolve();

function test(name, fn) {
  pending = pending.then(async () => {
    try {
      await fn();
      passed++;
      console.log(`  ok   ${name}`);
    } catch (err) {
      failed++;
      failures.push({ name, err });
      console.log(`  FAIL ${name}`);
      console.log(`       ${err.message}`);
    }
  });
}

function report() {
  return pending.then(() => {
    console.log(`\n${passed} passed, ${failed} failed`);
    if (failed > 0) process.exitCode = 1;
  });
}

module.exports = { test, report, assert };
