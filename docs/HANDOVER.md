# Handover — Clann MAG Explorer

Written for whoever (human or Claude) picks this up next. Read this first, then
[`docs/phase1-investigation.md`](phase1-investigation.md) for the detailed per-phase design log (every scope
decision, bug found, and performance number is recorded there — this document is the map, that one is the
territory). The original spec is [`clann-mag-explorer-brief.md`](../clann-mag-explorer-brief.md) at repo root.

## Where things stand

Phases 1–7 of the brief's 9-phase plan are done. All 126 tests pass (`node test/run.js`), the working tree is
clean, and the app has been verified end-to-end in a real browser at every phase (not just unit-tested).

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
| 8 | Comparison views and QC across the full set | **Not started** |
| 9 | Export and site chrome (About/FAQ, responsive pass) | **Not started** |

Run `node test/run.js` to confirm nothing regressed before you start. Every commit on `main` so far is one phase
(`git log --oneline` shows the full sequence) — keep that convention if it still makes sense.

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
  `app.js` set only inside `working-assignment.js`'s single mutation primitive — not on filter/search state
  changes (there are none yet; if Phase 8/9 adds any, decide deliberately whether they should also arm it, per
  the brief's own open question on this).
- **No coverage-aware outlier signal unless a coverage table was loaded**, and even then only when *every*
  contig in a bin has depth data — falls back to composition-only rather than partially scoring a mix. See
  `outliers.js`.

## What's next: Phase 8 (comparison views and QC across the full set)

Per the brief:
- Good-bin vs. bad-bin side-by-side view (tight, well-separated cluster vs. scattered/fragmented one) — likely
  reuses the Phase 7 scatter (`src/viz/scatter.js`, `scatter-geometry.js`) rendered twice with two different
  bins' contigs, since the plotting/hit-testing machinery already exists.
- Completeness vs. contamination scatter across *all* putative MAGs, coloured by which tool(s) support each one
  — another natural fit for the existing scatter code, this time one point per MAG rather than per contig
  (`bin-summary.js`'s `computeBinSummaries` / `bin-reconciliation.js`'s `putativeMags` already have everything
  needed: completeness, redundancy, and which tools' bins contributed).
- Redundancy check across finished MAGs to flag likely duplicate genomes split across separate bins — this is
  new logic, not yet built anywhere: probably a pairwise composition/marker-family-overlap similarity between
  MAGs, flagging pairs that look like the same organism binned twice.

None of Phase 8's three bullets are implemented yet. Skim the relevant section of `clann-mag-explorer-brief.md`
before starting, and check `phase1-investigation.md`'s existing findings sections for the established patterns
(pure-logic-first-then-DOM-glue, scope decisions recorded explicitly, real bugs caught by working through the
math or testing in-browser) before writing new code — the project has a consistent voice and testing discipline
by now and it's worth matching rather than reinventing.

## Then Phase 9 (export and site chrome)

- Revised contig→bin assignment table export (from `working-assignment.js`'s current state — `assignmentToRows`
  already produces the right shape, just needs a CSV/TSV serializer and a download trigger).
- Per-MAG FASTA extraction via `Blob.slice` against the Phase 2 `.fai`-style index (`fasta-index.js` already
  tracks byte offsets/line-layout per contig; the actual slice-and-reassemble-FASTA logic doesn't exist yet).
  Note the known limitation already documented in `phase1-investigation.md`: a gzip-compressed source's offsets
  are in the *decompressed* stream, so extraction from a `.gz` input needs a fresh decompression pass, not a
  direct slice — `faiEntry.sourceCompressed` is already there to detect this case.
- Summary table export per MAG (length, GC, N50, tool-agreement fraction, completeness tier) — mostly a CSV
  serializer over data `bin-summary.js`/`bin-reconciliation.js` already compute.
- About/FAQ content, responsive layout pass, and the `beforeunload`-adjacent "known limitation" framing the
  brief asks for (reload loses the live `File` reference — already a documented, accepted limitation per Phase
  1, not something to build around).

## A few things worth doing regardless of which phase you start with

- **The brief's open points section** (bottom of `clann-mag-explorer-brief.md`) lists several items marked as
  "a starting point... worth revisiting once tested against real student assemblies" — the paralog-safety
  thresholds, the pre-clustering aggressiveness, the chimerism-flagging rank cutoff. None of this has been
  tested against a *real* student assembly yet, only synthetic examples and the reference set's own held-out
  data. If real assemblies become available, that calibration work is still open.
- **No coverage-table or Kraken2-input real-world file has been tested**, only hand-built examples
  (`examples/marker-gene-demo.coverage.tsv`, `examples/marker-gene-demo.kraken2.tsv`). The MetaBAT2
  `jgi_summarize_bam_contig_depths` format handling in `coverage-table.js` was written from documentation, not
  verified against real tool output — worth a sanity check if a real file becomes available.
- **`git log` is your friend here** — every phase's commit message is a mini design doc in itself (rationale,
  what was tested, what was deferred). Read the last 2–3 before starting new work in an area you haven't
  touched yet.
