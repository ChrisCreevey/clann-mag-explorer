# Handover — Clann MAG Explorer

Written for whoever (human or Claude) picks this up next. Read this first, then
[`docs/phase1-investigation.md`](phase1-investigation.md) for the detailed per-phase design log (every scope
decision, bug found, and performance number is recorded there — this document is the map, that one is the
territory). The original spec is [`clann-mag-explorer-brief.md`](../clann-mag-explorer-brief.md) at repo root.

## Where things stand

All 9 phases of the brief's plan are done. All 145 tests pass (`node test/run.js`), the working tree is clean,
and the app has been verified end-to-end in a real browser at every phase (not just unit-tested).

| Phase | What it added | Status |
|---|---|---|
| 1 | Investigation doc, reference-set check, eDNA Explorer reuse audit | Done |
| 2 | Streaming FASTA parse, per-contig stats, `.fai`-style index | Done |
| — | Performance follow-ups (Web Worker, typed-array translation) | Done, ~11x speedup measured |
| 3 | Marker-gene search (seed-and-extend, BLOSUM62, paralog safety) | Done |
| — | Per-family threshold calibration | Done |
| 4 | Single contig→bin table loading, per-bin summaries, MIMAG tier | Done |
| 5 | Cross-tool bin matching/reconciliation, core/disputed sets | Done |
| 6 | Outlier flagging (composition/coverage/marker/taxonomy/Kraken2) | Done |
| 7 | Interactive reassignment (scatter + drag-select + live stats) | Done |
| 8 | Comparison views and QC across the full set | Done |
| 9 | Export and site chrome (About/FAQ, responsive pass) | Done |

There is no more scoped phase work from the brief. If you're picking this up next, the likely next step is either
(a) the calibration/real-data follow-ups in "Open items still worth doing" below, or (b) whatever new request
brought you here — run `node test/run.js` first either way to confirm nothing regressed.

## Architecture, in one paragraph

Vanilla HTML/CSS/JS, no build step, no runtime dependencies, GPL-2.0 — matches the other four Clann tools.
Every `src/parsers/*.js` and `src/model/*.js` module is a self-invoking function that attaches its exports to
`self` (not `window`), so the *identical* file loads via a plain `<script>` tag on the main thread and via
`importScripts()` inside `src/workers/fasta-worker.js` (`self === window` on the main thread, so one export
convention works both places). `src/app.js` is the only place DOM-specific code lives. `build/` is a separate
Node-based offline pipeline (allowed to use Node built-ins freely — it's never shipped) that turns
`reference-data/scg40_raw.fasta` into the static assets in `data/`.

## The five things most likely to bite you

1. **`<script>` tag order in `index.html` is a real dependency graph, not just a list.** Several modules
   destructure another module's exports at the top of their IIFE (`const { X } = self.ClannMAG.Y`), which means
   `Y`'s script tag must appear *before* the one that needs it. This bit me once already (`sniff.js` loaded
   before `breport.js`/`contig-bin-table.js`, which it depends on) — it passed every unit test (Node's
   `require()` graph doesn't care about order) and only broke in the actual browser. If you add a new module
   with a cross-module dependency, check `index.html`'s script order, and re-verify in-browser, not just via
   `node test/run.js`.

2. **The harness's `report()` prints a *cumulative* running total, not a per-file count.** `test/harness.js`'s
   `passed`/`failed` counters are shared across every file `test/run.js` `require()`s (same module instance,
   cached by Node). So when you see "48 passed" after one file and "69 passed" after the next, that's not a bug
   — it's 48 total so far, then 21 more. Don't be alarmed by numbers that look too high per file; check the
   final total.

3. **When testing in the Claude Browser tool against a local `python3 -m http.server`, the browser can serve a
   stale cached copy of a JS file even after a fresh `navigate()`, even in a brand-new tab on the same
   profile.** If you edit a file and the browser doesn't seem to see the change (e.g. an old error message keeps
   showing, or `window.ClannMAG.someExport.toString()` still shows old code), don't trust it — kill the Python
   server and restart it on a **different port**. That reliably busts the cache. Cost me a confusing debugging
   detour in the Phase 5 session; don't repeat it.

4. **Verify the actual math before trusting a test that passes.** Two real bugs this project shipped and then
   caught were both *logic* bugs, not typos, and both had tests that were passing right up until I worked
   through the arithmetic by hand:
   - Phase 5: a hand-built 3-tool "disagreement" test scenario turned out not to test disagreement at all,
     because a bin with only one competitor is trivially its own reciprocal best match regardless of how weak
     the overlap is — the test's premise was wrong, not the code (though the code's real behavior was also worth
     understanding). Fixed by recomputing the scenario so the dissenting bin's Jaccard overlap genuinely failed
     the threshold.
   - Phase 6: the first per-contig marker-taxonomy "distance from consensus" metric combined a contig's own
     taxIDs with the *whole-bin* consensus LCA and measured the gap. That's provably always zero — the whole-bin
     LCA is by construction an ancestor of every subset's LCA, so combining it back in in is a no-op. I verified
     this with a 3-line Node script (`t.lca([11, consensus])` really does just return `consensus`) *before*
     writing the "fixed" version, rather than trusting a test I hadn't yet written to catch it. Now uses a
     leave-one-out design (compare against the *other* contigs' consensus, not the whole bin's).
   Moral: when a per-contig or per-pair metric involves combining a subset with an aggregate of that same
   subset, stop and ask whether the operation can be a no-op by construction before writing the test.

5. **A stray literal `#` where a `//` comment was intended will break the file, and it's happened more than
   once** (auto-generated by me mid-edit, not a paste error) — e.g. a `#` at the start of a continuation line in
   a multi-line `//` comment block. `node -c <file>.js` catches it instantly; I now run that after any
   multi-line comment edit before moving on. If a file mysteriously fails to parse after an edit, check for this
   first — `grep -n "^#" src/**/*.js` finds it fast.

## Scope decisions you should know about (all deliberately made, all documented in `phase1-investigation.md`)

- **`build/03-calibrate.js`'s per-family thresholds**: only 1 of 40 marker-gene families (COG0495) needed an
  override. Don't assume the rest need tuning too — this was checked against all 40, not assumed from a
  two-family spot check (an earlier spot check *did* wrongly suggest COG0016 was fragile; it wasn't, once tested
  properly against 150 held-out sequences per family).
- **Kraken2 input parser is *not* the eDNA Explorer's `.breport` parser**, despite the brief's literal wording
  ("reuses the eDNA Explorer's Kraken2 report parser"). `.breport` is a whole-run aggregate with no per-query
  field — it structurally cannot answer "what did Kraken2 call this one contig." `src/parsers/kraken2-contigs.js`
  parses Kraken2's *other* standard output (`--output`, one row per query) instead. The disagreement check built
  on it is exact-taxID majority-vote, not lineage-aware, because a lineage-aware version would need a taxonomy
  tree covering whatever arbitrary taxa Kraken2's full database can call — out of scope, unlike the marker-gene
  provenance check, which only ever touches the ~5,000 taxa the fixed 40-family reference set references.
- **Interactive reassignment uses a rectangular drag-select, not a freeform lasso.** Recorded explicitly in
  `scatter-geometry.js`'s header. A real lasso is point-in-polygon hit testing — meaningfully more work for the
  same "select a visual cluster" outcome in a teaching tool.
- **`data/scg40-lineage.json` only covers taxa the marker-gene reference set itself touches** (1,753 referenced
  taxIDs + ancestors = 5,028 nodes, resolved from a real NCBI taxdump download in `build/03-taxonomy.js`), not
  the full NCBI taxonomy. It's fine for marker-gene provenance lineage; it is *not* usable for general Kraken2
  taxID lookups (see previous point).
- **The `beforeunload` guard arms on the first *reassignment* specifically**, via one module-level flag in
  `app.js` set only inside `working-assignment.js`'s single mutation primitive — not on filter/search state.
  Phases 8/9 added no new mutable session state that changes what's on disk if lost, so this was kept as-is.
- **No coverage-aware outlier signal unless a coverage table was loaded**, and even then only when *every*
  contig in a bin has depth data — falls back to composition-only rather than partially scoring a mix. See
  `outliers.js`.
- **Phase 8's MAG-redundancy check uses composition similarity, not marker-family overlap**, to flag likely
  duplicate genomes split across bins. The 40-family reference set is universal single-copy genes, so two
  unrelated MAGs would both hit most of the same families regardless of species — family identity carries no
  organism-identity signal once you're comparing *across* bins rather than completeness *within* one. Composition
  (mean tetranucleotide vector + GC) is the same signal `outliers.js` already uses to separate organisms within a
  bin; reused here in `src/model/mag-redundancy.js` for consistency. See `phase1-investigation.md`'s Phase 8
  section for the full reasoning and the in-browser verification (a synthetic near-identical-composition MAG
  pair correctly flagged, a compositionally distinct third MAG correctly not flagged).
- **Phase 9's per-MAG FASTA export skips contigs whose FASTA lines weren't uniformly wrapped**
  (`faiEntry.uniform === false`) rather than guessing a byte span for them — `src/model/fasta-extract.js`'s
  `computeSequenceByteSpan` returns `null` in that case and the caller surfaces a "N contig(s) skipped" note. A
  gzip-compressed source needs one full decompression pass before slicing (offsets are in the *decompressed*
  stream — `faiEntry.sourceCompressed`), done once per export via `DecompressionStream`, not per contig.
- **Export always reflects the live working assignment**, not the originally loaded tables — `initExportSection`
  in `app.js` reads straight from the same `workingAssignment` Map the Phase 7 interactive section mutates, via
  `assignmentToRows`. This matches the brief's "revised... table, reflecting any manual reassignment" wording
  literally: there's deliberately no separate "export the original" option.

## Open items still worth doing (not scoped phase work — pick up only if asked)

- **The brief's open points section** (bottom of `clann-mag-explorer-brief.md`) lists several items marked as
  "a starting point... worth revisiting once tested against real student assemblies" — the paralog-safety
  thresholds, the pre-clustering aggressiveness, the chimerism-flagging rank cutoff. None of this has been
  tested against a *real* student assembly yet, only synthetic examples and the reference set's own held-out
  data.
- **No coverage-table or Kraken2-input real-world file has been tested**, only hand-built examples
  (`examples/marker-gene-demo.coverage.tsv`, `examples/marker-gene-demo.kraken2.tsv`). The MetaBAT2
  `jgi_summarize_bam_contig_depths` format handling in `coverage-table.js` was written from documentation, not
  verified against real tool output.
- **The Phase 9 FASTA export hasn't been tested against a gzip-compressed source** — the decompression-then-slice
  path (`extractBinFasta`'s `needsDecompression` branch in `app.js`) is implemented per the documented offset
  limitation but only exercised by the uncompressed path so far, since the shipped example FASTA isn't gzipped.
- **No zip/multi-file bundling for FASTA export** — Phase 9 exports one bin's FASTA per click (a `<select>` +
  button), not "export all bins at once." Deliberately kept simple (no new dependency, matches the "no runtime
  dependencies" convention) rather than out of any brief requirement; revisit only if a real workflow needs it.
- **`git log` is your friend here** — every phase's commit message is a mini design doc in itself (rationale,
  what was tested, what was deferred). Read the last 2–3 before starting new work in an area you haven't
  touched yet.
