# Phase 1 investigation — data model and build pipeline plan

Status: draft for review, no code written against this yet except the reference-set check below.

## 1. Reference set — confirmed against the brief

`reference-data/scg40_raw.fasta` (formerly `assets/all.fas`) matches the brief exactly:
- 69,208 sequences
- 40 distinct COG families
- Header format `>COGID.fa.<taxID>.<locus_tag>`, e.g. `>COG0012.fa.1000565.METUNv1_03812`

## 2. eDNA Explorer reuse — what's actually there

Cloned `chriscreevey/clann-edna-explorer` for reference. Findings that change what "port the LCA logic" means in practice:

- **Architecture matches the brief's expectation**: no build step, per-module plain `<script>` tags (`src/parsers/*.js`, `src/model/*.js`, `src/viz/*.js`), each module attaches to `window.ClannEDNA.*`, also exports via `module.exports` for its own test harness (`test/*.test.js` + `test/run.js`, no external test framework). We should follow this exact pattern.
- **`src/model/taxonomy-tree.js`** has a reusable `TaxonomyTree` class: flat, typed-array-backed (`taxid[]`, `parentIndex[]`, `depth[]`, `rankLetter[]`/`rankSub[]`), taxid→dense-index map, `getOrCreateNode`, per-node `perSample` counts. This data structure is directly reusable for a taxID→lineage table.
- **No LCA function exists anywhere in that codebase.** The brief's phrasing ("port the tree-traversal/LCA logic") implies a function to lift; there isn't one — only the tree data model. What actually needs writing fresh: (a) an LCA-of-N-taxa function walking `parentIndex` to root and finding the deepest common ancestor, and (b) a "rank at which lineages agree" helper for the chimerism flag. Both are small (~30 lines) but are new code, not ported code.
- **The tree is currently populated from a Kraken breport's indentation structure** (parent inferred from the tree's local depth as you read report lines in order). Our tree instead needs to be built from `nodes.dmp` (taxid, parent taxid, rank) + `names.dmp` (taxid, scientific name) in an NCBI taxdump — a different, and simpler, ingestion path (flat parent-pointer table, no indentation to interpret).
- **Kraken2 report parser** (`src/parsers/breport.js`) is genuinely reusable as-is for the optional per-contig taxonomic-call input, per the brief.

**Conclusion**: reuse the `TaxonomyTree` class and the breport parser verbatim (adapted for per-contig rather than per-sample counts); write the taxdump ingestion and LCA logic fresh, informed by the tree's existing shape so the two stay compatible.

## 3. Streaming FASTA + index strategy

Plan, to be validated against a real multi-hundred-MB assembly before Phase 2 starts:

- Read via `file.stream().getReader()` in chunks (e.g. 1–4 MB), decode as UTF-8 incrementally (a `TextDecoder` with `{stream: true}`), buffer only the current partial line across chunk boundaries.
- Track byte offset as we go. On each `>` header line: close out the previous record (store its `.fai`-style entry — offset of first sequence byte, sequence length in bases, bases-per-line, bytes-per-line, assuming uniform line wrapping as `samtools faidx` does — and fall back to per-line bookkeeping if a contig's wrapping is irregular), then start accumulating the new record's stats.
- Per-contig running stats computed incrementally per chunk (no need to hold the full sequence): length, GC count, GC-skew running sums, tetranucleotide k-mer counts (a `Uint32Array` or plain object keyed by 4-mer, ~136 canonical counts after reverse-complement collapsing), six-frame translation done incrementally is awkward across chunk boundaries for reading frames 2/3 — simplest correct approach is to buffer the current contig's raw sequence only (not the whole file) until its terminating `>` or EOF, then translate/search it, then discard. This is still "one contig in memory at a time," consistent with the brief's constraint (assemblies have short-to-moderate contig lengths compared to the whole assembly; buffering one contig, even a large one, is orders of magnitude cheaper than buffering the whole file).
- `File`/`Blob` reference retained for the session; sequence extraction later uses `blob.slice(start, end)` + decode, driven by the `.fai`-equivalent index.
- Decision carried over from your answer above: **no persistence across reload** — accept the limitation, document it in the UI, no File System Access API handle work in Phase 1.

Open sub-question I'd like to confirm: tetranucleotide (4-mer) vs. a different k for the composition signature — 4-mer is the standard choice in binning literature (matches MetaBAT2/CONCOCT-style composition vectors) and I'd default to it unless you want something else.

## 4. Marker-gene offline build pipeline

Three steps, in order, each a separate Node script under `build/`:

1. **`build/01-cluster.js`** — pre-cluster `reference-data/scg40_raw.fasta` within each COG family at ~90% identity (single-linkage greedy clustering via a fast approximate method — minhash or straightforward k-mer-profile similarity, since exact all-pairs alignment across ~1,700 sequences/family × 40 families is unnecessary for a 90% cutoff). Output: `build/intermediate/scg40_clustered.fasta` + a representative-count-per-family log for eyeballing how much diversity survives.
2. **`build/02-index.js`** — build the seed index from the clustered set: reduced-alphabet translation (Murphy-10 or similar 10–12 group scheme — I'll pick a published reduced alphabet rather than invent one), k=5 windowing (k=6 as a fallback if k=5 produces too many index collisions during calibration), hash to integer keys, build parallel typed arrays (`Int32Array` key table + offset/value arrays, CSR-style) rather than a JS `Map`. Output: `docs-site/data/scg40-index.bin` (binary) + `docs-site/data/scg40-refseqs.bin` (real residues for extension scoring, also typed-array/binary, not FASTA text).
3. **`build/03-taxonomy.js`** — download current NCBI taxdump (`taxdump.tar.gz`) at build time only (not checked into the repo), extract the distinct taxIDs actually present in `scg40_raw.fasta` headers, resolve each against `nodes.dmp`/`names.dmp`, handle merged/deleted taxIDs via `merged.dmp`/`delnodes.dmp` (map merged IDs forward, drop or flag dead ones with a warning), output a compact `docs-site/data/scg40-lineage.bin` (or JSON if small enough — few thousand taxa, likely fine as JSON, will confirm size during calibration) covering only the referenced taxa.

All three scripts are checked into the repo (`build/`), runnable via plain `node build/0N-*.js`, no bundler. Their outputs (the `.bin` assets) are also checked into the repo per your earlier framing in the brief ("versioned so they can be regenerated"), so the deployed site never needs the build step to run.

## 5. Binary index format (concrete proposal)

- `scg40-index.bin`: header (magic bytes, format version, k, alphabet-size) + CSR-style arrays: `keyOffsets: Uint32Array` (one entry per possible reduced-alphabet key, or a sorted-key + binary-search scheme if the key space is sparse — will measure during calibration which is smaller), `hitRefSeqId: Uint16Array`, `hitPosition: Uint16Array`.
- `scg40-refseqs.bin`: header + `seqOffsets: Uint32Array` + concatenated residue bytes (`Uint8Array`, one byte per amino acid, real alphabet) + a small side table of (refSeqId → COG family index, taxID).
- Loaded via `fetch().arrayBuffer()` + typed array views, no parsing step beyond reading the header.

This will get refined once real numbers exist (index size, collision rate) — flagging it as a plan, not a commitment, until calibration runs.

## 6. Threshold calibration plan

Deferred to actual implementation once the pipeline exists (per your answer: build first is being sequenced as design → code, so this is downstream of steps 1–5 above, not blocking them). Will test seed/extend thresholds against six-frame translations of a genome (or a few) not expected to be enriched for these 40 families, per the brief.

## Phase 2 findings — streaming FASTA parsing and per-contig stats (implemented)

Built per the plan above: [src/model/fasta-index.js](../src/model/fasta-index.js) (streaming reader +
`.fai`-style index, gzip-aware), [src/model/translate.js](../src/model/translate.js) (six-frame translation,
sentinel-delimited, shared with the future marker-gene module per the brief), and
[src/model/contig-stats.js](../src/model/contig-stats.js) (GC/GC-skew/tetranucleotide composition/coding-density).
26 unit tests pass (`node test/run.js`), and the load path is wired into the UI end to end (verified in-browser:
loading a file renders the assembly summary + per-contig table with no console errors).

**Performance finding, flagged rather than silently shipped**: streaming a realistic synthetic 50 MB assembly
(4,739 contigs, 500bp–20kb each) end-to-end — parsing, `.fai` index, GC/skew/composition, and six-frame
translation for the coding-density estimate — took **22.5s** (~2.2 MB/s), dominated by the six-frame translation
step (confirmed via isolated benchmarking: translation is roughly 3-5x the cost of everything else combined).
One round of optimization already applied (array-join instead of repeated string `+=` in
`reverseComplement`/`translateFrame`, cutting translation time ~43%); further gains would need a deeper rewrite
(numeric/typed-array codon encoding instead of JS string ops) that I haven't done, given this is exactly the
computation the brief already earmarks as reused, at the same cost, by Phase 3's marker-gene search — so this
number is a preview of Phase 3's cost too, not just Phase 2's.

**Consequence for the "hundreds of MB" assemblies the brief anticipates**: at this throughput, a 300 MB assembly
would take on the order of two minutes, blocking the main thread the whole time (no progress UI, no
responsiveness) since nothing here runs off-thread yet. Worth a decision before Phase 3 compounds the cost:
- Move the streaming parse + translation to a Web Worker, so the page stays responsive and a real progress bar
  is possible — doesn't reduce total time, but changes "frozen tab for 2 minutes" into "visible progress".
- And/or optimize `translateFrame`/`reverseComplement` further with typed arrays and numeric codon lookup instead
  of string concatenation and object-key hashing — a genuine speed win, not just perceived responsiveness.
- Both are worth doing before Phase 3 (marker-gene search) adds the seed-and-extend cost on top of the same
  translation.

Not blocking further Phase 2/3 development, but flagging now rather than after Phase 3 makes it worse.

**Update — both follow-ups implemented before starting Phase 3:**
- **Typed-array translation** ([src/model/translate.js](../src/model/translate.js),
  [src/model/contig-stats.js](../src/model/contig-stats.js), new shared
  [src/model/dna-codes.js](../src/model/dna-codes.js)): per-base string allocation and object-key hashing
  replaced with 2-bit base codes, a 64-entry numeric codon table, a rolling integer k-mer code for the
  composition scan, and a single `TextDecoder` pass instead of per-character string building. Public API
  unchanged (string in, string out), so no test or caller changes needed beyond re-running them.
  **Result: the same 50MB/4,773-contig benchmark dropped from 22.5s to 3.17s — a 7x speedup.** A 300MB assembly
  is now on the order of ~19s of compute instead of ~2 minutes.
- **Web Worker** ([src/workers/fasta-worker.js](../src/workers/fasta-worker.js)): the streaming parse now runs
  off the main thread. `app.js` posts the `File` to the worker, which streams it via the same `fasta-index.js`
  and posts per-contig records plus progress pings (every 200 contigs) back; the UI renders a live "N contigs so
  far" message and stays responsive throughout, instead of freezing for the parse's duration. Verified in-browser
  with a 900-contig/3.1MB synthetic file (0.2s) and the smaller example assembly — no console errors either way.
- Both changes required switching every model/parser module's environment guard from `window` to `self`
  (`self === window` on the main thread, but only `self` exists inside a Worker), so the same unmodified files
  load via both `<script>` tags and `importScripts()`. This also means Phase 3's marker-gene search can reuse
  these modules directly inside a worker (the same or a new one) without further changes.

**Second typed-array follow-up (user-suggested):** the reverse-frame translation was still building a full
reverse-complement string, then re-scanning it three times with the same base-code lookup already paid for once
during that construction — a genuine redundancy, since a codon's reverse-complement amino acid is a fixed
function of that codon alone. Replaced with a 64-entry `RC_CODON_CHAR_CODE` table (forward codon code →
reverse-complement translation) and `translateReverseFrame(seq, offset)`, which walks the *forward* sequence
right-to-left in triplets and does one table lookup per codon — no reverse-complement string ever built.
`reverseComplement()` itself stays exported (still useful for Phase 9 sequence extraction), just no longer used
inside `translateSixFrames`.

Verified equivalent to the old `translateFrame(reverseComplement(seq), offset)` behaviour via a property test
across sequence lengths 0–1000 (including lengths shorter than one codon and lengths not divisible by 3) and all
three offsets — this also caught a latent bug in both `translateFrame` and the new function: `Math.floor` on a
negative `(length - offset)` produced a negative typed-array length for very short sequences, now clamped to 0.

**Result**: translation alone is 4-9x faster in isolation (JIT-warmup-dependent), but the full pipeline only
improved 3.17s → 2.95s (~7%) on the same 50MB benchmark. This turned out to be premature — see below, actual
profiling showed translation was *still* the dominant cost.

**Third round — unified base-code pipeline (user-suggested, following the same reasoning as the RC-codon
table).** Rather than guessing further, profiled `computeContigStats` directly (50MB/~4,800-contig benchmark,
2.69s of the pipeline's 2.95s total):

| Step | Time | Note |
|---|---|---|
| `toUpperCase()` | 47ms | cheap alone, but unconditional every contig |
| GC-counting (`ch === 'G'` string comparison) | 520ms | surprisingly expensive vs. integer ops |
| Composition scan | 382ms | already typed-array based |
| Six-frame translation | 1,677ms | **62% of the total — still dominant**, correcting the "no longer dominant" guess above |

The real inefficiency: GC-counting, composition, and all six translation frames were each independently
re-deriving "which base is this" from the string — via three different techniques (string equality, `charCodeAt`
+lookup, `charCodeAt`+lookup again per frame) — touching most bases 6+ times combined.

**Fix**: [src/model/dna-codes.js](../src/model/dna-codes.js) gained `computeBaseCodes(seq)`, converting the whole
sequence to a 2-bit-per-base `Int8Array` once, with `BASE_CODE` extended to cover lowercase too (removing the
`toUpperCase()` pass entirely — no more separate case-folding copy). GC-counting and the composition scan are now
fused into one pass over that array in [contig-stats.js](../src/model/contig-stats.js)'s new `scanBaseCodes()`.
[translate.js](../src/model/translate.js) gained `*FromCodes` variants (`translateFrameCodes`,
`translateReverseFrameCodes`, `translateSixFramesFromCodes`) that read codons via pure integer array indexing
instead of `charCodeAt`+lookup; the original string-in/string-out functions (`translateFrame`,
`translateReverseFrame`, `translateSixFrames`) became thin wrappers for callers that don't already have a code
array, so the public API and all prior tests were unaffected.

Verified via existing tests plus two new ones (lowercase-vs-uppercase equivalence in both `contig-stats.test.js`
and `translate.test.js`, since the case-insensitive `BASE_CODE` is new behavior worth locking in) — 31 tests
pass. Prototyped the approach with real numbers before implementing: `computeContigStats` dropped from 2.69s to
1.52s (~44% further cut) in the prototype, and after full implementation the **complete `streamFasta` pipeline
dropped from 2.95s to 2.09s (~29%)** on the same 50MB benchmark.

**Running total across all three optimization rounds: 22.5s → 2.09s, a ~10.8x overall speedup**, verified
correct at each step (unit tests + in-browser check against the example assembly, same output every time) rather
than trading correctness for speed.

**Considered and rejected**: switching `current.seq += line` (accumulating a contig's sequence while streaming)
to array-push-then-join, on the general JS-performance folklore that repeated `+=` is slow. Measured it directly
instead of assuming — `+=` was actually faster on a representative sample (V8's rope-based concatenation handles
this pattern well) — so left unchanged. Not guessing where a quick measurement settles it.

**Still on the table, not done**: a byte-level rewrite of the line-parsing loop itself (operating on raw
`Uint8Array` chunks instead of decoding to strings, only decoding headers). Line parsing is a small slice of
total time (~9% before this round, likely less now that translation shrank), so this is a smaller/riskier-return
change than what's been done — worth it only if more throughput is needed later.

## Still open / needs your call before I start writing code

1. **Reduced alphabet scheme** — I'll default to a published one (e.g. Murphy10) unless you have a preference.
2. **Clustering method for step 1** — greedy approximate (minhash/k-mer-profile) vs. exact — I'd default to approximate for speed, willing to revisit if it under- or over-clusters during calibration.
3. **taxdump merged/deleted taxID handling** — default plan above (map merged forward, warn-and-drop truly dead ones) — confirm this is acceptable rather than something that should hard-fail the build.
4. Ready to proceed to repo scaffolding (site skeleton + `build/` stubs) once you've reviewed this, unless you'd rather adjust something above first.

## Phase 3 findings — marker-gene identification module (implemented)

Scope decision made and flagged rather than silently assumed: the brief's Phase 3 description ("produce per-contig
gene-family tags and provenance taxIDs") only needs each hit's *raw* taxID, not the full lineage table — the
taxID→lineage table and LCA/consensus-lineage computation are Phase 6 concerns (bin-level chimerism check) layered
on top of these per-contig tags once bins exist. So **`build/03-taxonomy.js` is deferred until Phase 6 actually
needs it**, not implemented as part of this phase — avoids building an asset nothing consumes yet.

**Pipeline built and run against the real reference set**, not stubbed:
- [build/01-cluster.js](../build/01-cluster.js): MinHash-based greedy clustering at the ~90%-identity-proxy
  threshold from the earlier plan. Real run: 69,208 → 39,854 representatives (57.6% retained), 15s.
- [build/02-index.js](../build/02-index.js): reduced-alphabet (Murphy10) seed index + real reference sequences.
  Real run output: **18.9MB index + 13.0MB reference sequences ≈ 32MB total shipped assets.**
- [src/model/reduced-alphabet.js](../src/model/reduced-alphabet.js), [blosum62.js](../src/model/blosum62.js),
  [marker-genes.js](../src/model/marker-genes.js): shared Murphy10 windowing (build and runtime must window
  identically or seeding can't find matches), the standard BLOSUM62 matrix, and the seed-and-extend search
  itself — seeding, diagonal binning, ungapped X-drop extension, and all three paralog-safety checks from the
  brief (score/coverage threshold, family-margin, multi-representative agreement) applied together.
- Wired end-to-end: [src/workers/fasta-worker.js](../src/workers/fasta-worker.js) loads the ~32MB of assets as
  soon as the worker starts (overlapping network time with the user picking a file), and runs marker search
  against each contig's six-frame translation — the same translation `contig-stats.js` already computes for
  coding density (`frames` now part of its return value), not recomputed. `app.js` shows a per-contig "Marker
  genes" column plus assembly-level "contigs with marker genes" / "distinct families hit" summary stats.

**Two real bugs caught by testing against actual data, not by the unit tests in isolation:**
1. Reference headers are `COGID.fa.<taxID>.<locus_tag>` — `parts[1]` is the literal string `"fa"`, not the
   taxID (`parts[2]` is). The initial `parseHeader` read `parts[1]`, silently producing `taxId=0` for every one
   of the 39,854 representatives. Unit tests using hand-built fixtures didn't catch this (the fixtures had
   correct taxIDs baked in already); running the real built index against a real contig and checking the
   reported provenance taxID did. Fixed, and a regression test now covers `parseHeader` directly
   ([test/build-index.test.js](../test/build-index.test.js)).
2. A latent typed-array-length bug carried over from Phase 2's translation work (see above) — already fixed
   there, re-verified it doesn't recur here since `marker-genes.js` reuses `translateSixFramesFromCodes` rather
   than reimplementing translation.

**Performance finding, the most consequential one in this phase — measured, not assumed:** an early version
(k=5 reduced-alphabet seeding, no per-diagonal seed-count gate, 100-hits-per-key cap) made *every* window of
even a fully random, marker-free contig a candidate for expensive extension: 81.5% of the 100,000 possible k=5
reduced-alphabet keys are populated in this reference set (real protein sequences aren't remotely uniform over a
10-letter alphabet), so a 20kb random contig took **1.6 seconds** to search — projecting to roughly an hour for
a realistic 50MB assembly. Three fixes, each measured before/after rather than guessed:
- **k=6 instead of 5** (1,000,000 possible keys): cuts populated-key occupancy to ~45%. 20kb random: 1.6s → 0.58s.
- **Two-hit heuristic** (BLAST's classic fix for exactly this problem): only extend a `(refSeqId, diagonal)`
  once it accumulates ≥2 distinct seeds, not on the first isolated one — a real homologous region produces many
  seeds on the same diagonal; one isolated seed is usually chance noise. 20kb random: 0.58s → 0.32s.
- **Lower per-key hit cap** (100 → 20; median is only 5 hits/key at k=6, so this leaves the large majority of
  keys completely untouched): directly cuts the seeding/bookkeeping cost specifically for the generic,
  over-represented keys a random sequence is most likely to hit. 50kb random: 0.80s → **0.24s**.

**Net result: roughly 16x faster than the naive version, real detection unaffected** (verified against real
COG0012/COG0016 fragments embedded in synthetic contigs throughout — both still call correctly with the tighter
parameters). **Honestly, still not fast**: 0.24s per 50kb extrapolates to roughly **4 minutes for a realistic
50MB assembly** — spent inside the same Worker that already keeps the page responsive during parsing (Phase 2),
so this is a "the load takes several minutes, with visible per-contig progress" problem, not a frozen-tab
problem, but still a real UX cost worth revisiting rather than quietly accepting. **Not done, flagged as
follow-up** rather than pursued further given the scope already spent on this phase:
- A minimizer-style seed selection (only check e.g. every 3rd window as a sketch, instead of every position) —
  a genuine algorithmic speedup, not just parameter tuning, but a bigger change to get right.
- More aggressive clustering (lower the ~90%-identity threshold further) — directly shrinks the index and
  therefore seeding cost, at a real diversity cost to the multi-representative-agreement check.
- Calibrating `MIN_SEEDS_PER_DIAGONAL`, `MAX_HITS_PER_KEY`, and `k` together against real assemblies (the
  brief's own Phase 1 calibration ask) rather than the random-sequence proxy used here.

**Calibration placeholders, explicitly not validated** (in `marker-genes.js`'s `DEFAULT_PARAMS`): `xDrop=15`,
`minScore=50`, `minCoverage=0.5`, `minMargin=10`, `minRepresentatives=2`. These are reasonable starting points,
not numbers tested against negative-control genomes as the brief's Phase 1 calibration explicitly calls for —
exposed as overridable `params` specifically so real calibration doesn't require code changes later.

**Tests**: 12 new tests across [test/marker-genes.test.js](../test/marker-genes.test.js) (binary round-trip,
extension scoring and sentinel-halting, real-sequence family discrimination between two actual different COG
families, negative control, both directions of the multi-representative-agreement check) and
[test/build-index.test.js](../test/build-index.test.js) (the taxID-parsing regression). 40 tests total, all
passing. Verified in-browser against a synthetic assembly with two contigs each containing a real, different
embedded marker gene plus one marker-free contig — all three called correctly, no console errors.

## Phase 3 calibration findings — per-family thresholds, not one global guess

Prompted by a direct question worth recording: a hand-picked two-family spot check (COG0012 vs COG0016, used
throughout the performance-tuning work above) showed wildly different natural headroom — COG0012 tolerated
`MIN_SEEDS_PER_DIAGONAL` up to 6 with 95 representatives to spare, COG0016 lost its call entirely past 3. That
was real evidence the *single global* `DEFAULT_PARAMS` (guessed placeholders per the note above) could be
simultaneously too strict for a thin family and needlessly loose for a robust one — but two hand-picked families
is not "checked against all 40," and the fix belongs in calibration data, not another guess.

**The held-out set already existed for free.** `build/01-cluster.js` discards ~42% of the raw reference set as
near-duplicates of a chosen representative — real, known-family sequences that never went into the shipped
index. Extended it to also write `build/intermediate/scg40_heldout.fasta` (29,354 sequences) instead of
discarding them silently: a natural true-positive test set with no external genome download needed.

**[build/03-calibrate.js](../build/03-calibrate.js)**: for each of the 40 families, samples 150 held-out
sequences (evenly strided, deterministic), runs each through the *same* search code path production uses
(`computeFamilyCandidates` — a refactor of `searchContigForMarkers` that exposes the pre-filter candidate data:
representative count and margin-over-next-best-family *before* the margin/agreement gates are applied, so
calibration can see what a query could have supported, not just whether it happened to clear today's
threshold). Derives, per family, the 10th-percentile representative-count and margin observed across its own
true positives — i.e. the largest threshold ~90% of that family's real members would still clear — discounted
by 20% for safety margin beyond the sample. Deliberately **only overrides a family when the default would cost
it measured sensitivity**, never tightens beyond default even where a family has spare headroom: the held-out
set is drawn from the same collection the index itself came from, so it's a trustworthy *relative* signal
between families but not a substitute for the brief's own ask of testing against genomes not enriched for these
families — tightening further based on this data alone would be fitting to the reference set's internal
diversity, not real-world divergence.

**Real result, run against the full shipped index**: only **1 of 40 families** (COG0495) needed an override
(`minRepresentatives` 2→1; its held-out true positives' 10th-percentile representative count was 2, discounted
to 1). Every other family's margin comfortably cleared the default `minMargin=10` by two to three orders of
magnitude (observed 10th-percentile margins ranged ~390–4,300 across families) — meaning the earlier
COG0016-at-threshold-4 finding was specific to the diagonal-seed-count knob and the one hand-picked fixture
sequence used to test it, not a sign the family is generally fragile; tested against 150 real held-out COG0016
members, it turned out to be unremarkable. This is the concrete payoff of measuring against all 40 rather than
extrapolating from two: the actual problem was much smaller and more localized than the spot check suggested.
Also surfaced, not acted on (out of scope — would mean calibrating `minScore`/`minCoverage` per family, which
carries the same reference-internal-diversity overfitting risk flagged above): two families (COG0124: 1/150,
COG0525: 3/150) had a handful of held-out sequences fail even the base score/coverage threshold outright — a
low, plausibly-acceptable-under-the-brief's-"under-calling is fine" rate, but a real one worth knowing about.

**Runtime wiring**: `loadMarkerGeneAssets` now also fetches `data/scg40-thresholds.json` (the calibration
script's output — `{defaultParams, familyOverrides}`), falling back to `DEFAULT_PARAMS` for every family if the
fetch fails, consistent with the module's existing "optional, degrade gracefully" framing. `searchContigForMarkers`
resolves each family's effective parameters as `DEFAULT_PARAMS < assets.thresholds < explicit params argument`
(most specific wins), so a caller can still override anything at call time for testing without touching the
shipped asset. 8 new tests (3 for the override resolution/precedence behavior, 5 for `build/03-calibrate.js`'s
`evenSample`/`percentile` helpers) — 48 tests total, all passing.

## Phase 4 findings — single binning-result loading and bin summaries

**Scope**: brief's Phase 4 is deliberately single-tool (one contig->bin table); cross-tool reconciliation
(matching equivalent bins across multiple tools' tables) is Phase 5. Built:

- [src/parsers/contig-bin-table.js](../src/parsers/contig-bin-table.js): parses the brief's "two-column tab/CSV"
  input. Delimiter (tab vs comma) and header presence are both detected from content, not assumed — checked
  against the two real conventions in circulation: DAS_Tool's `Fasta_to_Contig2Bin.sh` output (tab-delimited,
  no header) and CONCOCT's `clustering_gt1000.csv` (comma-delimited, `contig_id,cluster_id` header row).
- [src/parsers/sniff.js](../src/parsers/sniff.js): now a real content-based sniffer, distinguishing `breport`
  (exactly 6 tab fields per line, reuses `parseBreportLine`) from `contig-bin-table` (uniform 2-column shape).
  Doesn't route to `coverage-table` yet since that parser is still a stub — sniffing a format this repo can't
  parse would be worse than not detecting it.
- [src/model/bin-summary.js](../src/model/bin-summary.js): per-bin `computeBinSummaries` joining loaded contig
  records to a bin table by ID — contig count, total length, N50/L50, mean GC, and SCG-based
  completeness/redundancy computed directly from Phase 3's per-contig `markerHits` (completeness = fraction of
  the 40 families found anywhere in the bin; redundancy = fraction of *found* families that appear on more than
  one contig, the standard CheckM-style extra-copy proxy for contamination), plus a MIMAG-style tier
  (`high`/`medium`/`low`) from completeness/contamination thresholds only — explicitly not the full MIMAG
  standard, which also needs rRNA/tRNA presence this module has no signal for; said so in the UI, not just here,
  per the marker-gene module's own precedent of surfacing limitations next to the number itself. Thresholds are
  an overridable object (`DEFAULT_MIMAG_THRESHOLDS`), matching the "shown and adjustable" framing in the brief.
  Contigs in the bin table with no matching loaded contig are reported (`unmatchedContigIds`), not silently
  dropped or thrown on — a mismatched ID is a plausible real mistake (wrong assembly/table pairing) worth
  surfacing.
- `app.js`: the existing multi-file picker now also content-sniffs every non-assembly file for a contig->bin
  table (first match wins — multiple tables is Phase 5) and renders a "Bin summaries" card between the assembly
  summary and per-contig table when one loads alongside the FASTA.

**One real bug caught by browser testing, not the unit tests**: `sniff.js` calls `parseBreportLine` from
`breport.js` and `looksLikeContigBinTable` from `contig-bin-table.js` at module-eval time (top-level
destructuring, not inside a function), so it needs both already loaded — `index.html` had it listed *first*
among the parser `<script>` tags, which threw immediately in the browser (`Cannot read properties of undefined
(reading 'breport')`) while every one of `sniff.js`'s own unit tests still passed, since Node's `require()`
graph doesn't care about the same ordering a browser's sequential `<script>` tags do. Fixed by reordering the
tags (`breport.js`/`contig-bin-table.js` before `sniff.js`); a reminder that this suite's "attach to `self`,
load via plain `<script>` tags" convention makes tag order a real correctness dependency the unit tests can't
see, the same category of gap Phase 3's `parseHeader` bug was, just at load time instead of parse time.

**Tests**: 21 new (`contig-bin-table.test.js`, `sniff.test.js`, `bin-summary.test.js`) — 69 total, all passing.
Verified in-browser: `examples/marker-gene-demo.bins.tsv` (added alongside the existing marker-gene demo FASTA)
loaded together correctly join, render the bin summary card, and show the right completeness numbers for the
two contigs' already-verified-real COG0012/COG0016 marker calls.

## Phase 5 findings — cross-tool bin matching and reconciliation

**The core feature** (brief §Background): works at contig granularity, not DAS_Tool's whole-bin granularity —
for each contig, how many of the loaded tools agree on where it belongs, computed by matching each tool's bins
to one another via contig overlap.

- [src/model/bin-reconciliation.js](../src/model/bin-reconciliation.js): bins from different tools are matched
  via **reciprocal-best-hit Jaccard overlap** — bin A (tool 1) and bin B (tool 2) are linked only if each is the
  other's single best match by contig-set Jaccard similarity AND that score clears `minJaccard` (default 0.1).
  Same conservative principle the marker-gene module's paralog-safety margin check uses, applied here so two
  only-loosely-overlapping bins don't get merged on weak evidence. Matched bins (across any number of tools, not
  just pairs) are grouped into putative MAGs via union-find over all pairwise reciprocal-match edges — makes
  matching transitive across 3+ tools without needing an explicit multi-way clustering step.
- **Per-contig agreement** comes from a simple vote: each tool that assigned a contig to *some* bin (not
  "unbinned" — see below) votes for whichever putative MAG that bin ended up in. `agreementFraction` =
  majority-vote-count / total-voting-tools. A contig is core to its majority MAG only at `agreementFraction ===
  1` (unanimous among tools that had an opinion); anything else goes to that MAG's `disputedContigIds` and,
  ranked by how split the vote is, `disputedContigsRanked`.
- **"Unbinned" handling, a real correctness issue caught before it became a bug, not after**: two different
  tools' "everything left over" buckets usually share most of their contigs by construction, which would make
  naive bin-matching merge them into one large bogus putative MAG. `UNBINNED_LABELS` (`unbinned`, `none`, `na`,
  etc., case-insensitive) excludes those rows from bin-matching entirely — the contig gets no vote from that
  tool, rather than a vote for a fake "everyone's leftovers" MAG. Deliberately conservative: only exact known
  labels, not e.g. any purely-numeric bin ID, so a tool's real bin "0" isn't swept in by mistake — documented as
  a known limitation, not exhaustive across every tool's own unbinned convention.
- **UI** (`app.js`): the file picker now content-sniffs *every* non-assembly file for a contig->bin table (not
  just the first, per Phase 4) and labels each by its filename (content alone can't name which tool produced a
  table). Two or more tables get a "Cross-tool reconciliation" card — a side-by-side putative-MAG table (each
  tool's contributing bin ID + contig count in its own column, core/disputed counts, and completeness/redundancy
  computed from the high-confidence core contig set only) plus a ranked disputed-contig table (each tool's vote
  shown per contig) — alongside the existing per-tool Phase 4 bin-summary cards, not instead of them.

**A real bug in my own test, not the implementation**, caught by actually running the numbers rather than
trusting the intent behind a hand-built scenario: a first attempt at a "3 tools, 2 agree vs 1 disagrees" test
put the disagreeing tool's alternative bin as a single-contig bin with nothing else to compare against — since
a bin with only one competitor is *trivially* its own reciprocal best match regardless of how weak the overlap
actually is, it got merged into the majority group by transitivity, silently defeating the disagreement the
test meant to construct. Fixed by padding the disagreeing bin with contigs no other tool mentions, driving its
Jaccard overlap with the majority group's bin decisively under `minJaccard` so it stays a separate MAG — a good
reminder that "give the minority voter its own bin" isn't sufficient to test disagreement; the bin's *overlap*
with the majority is what the algorithm actually keys on.

**Tests**: 9 new (`bin-reconciliation.test.js`) — 77 total, all passing. Verified in-browser with a second
example bin table (`examples/marker-gene-demo.bins2.tsv`, deliberately disagreeing with the existing one on
`contig_2`) loaded alongside the assembly and the first table: 2 putative MAGs matched correctly, `contig_2`
correctly flagged as the sole disputed contig at 50% agreement, both per-tool bin-summary cards and the
reconciliation card rendered with no console errors.

## Phase 6 findings — outlier flagging and taxonomic overlay

Five signals, per the brief, combined into one ranked per-contig view rather than five separate ones — matching
the brief's own framing ("a contig that is both a statistical outlier and cross-tool disputed is flagged as the
strongest reassignment candidate").

**`build/03-taxonomy.js`, deferred in Phase 3, implemented now that Phase 6 actually needs it**: downloads the
current NCBI taxdump (~78MB, build-time only), resolves the 1,753 distinct taxIDs referenced across
`reference-data/scg40_raw.fasta`'s headers plus every ancestor up to root (5,028 nodes total; 24 taxIDs had been
merged forward, 0 were dead), ships `data/scg40-lineage.json` (JSON, not binary — 419KB was well within the
"small enough" bound flagged as an open question back in Phase 1). `src/model/marker-taxonomy.js` builds a
`taxonomy-tree.js` `TaxonomyTree` from it (parent-before-child insertion order, since input rows aren't
topologically sorted) and reuses that class's `lca()` — the exact function Phase 3 wrote because the eDNA
Explorer's port had the tree shape but no LCA — to compute each bin's marker-provenance consensus lineage.

**A real math bug caught by working through the proof, not by a failing test**: the first version of the
per-contig taxonomic-distance metric combined a contig's own provenance taxIDs with the *whole-bin* consensus
taxID and measured how far that combination sat from the consensus. That's provably a no-op: the whole-bin LCA
is by construction an ancestor of every subset's LCA (including any single contig's), so adding it back into a
subset's own LCA can only ever return the consensus itself — verified directly (`t.lca([11, consensus])` really
does just return `consensus`) before it shipped as dead code with zero actual signal. Fixed with a **leave-one-
out** design instead: compare each contig's own provenance against the LCA of *every other contig in the bin*,
which isn't fixed relative to the contig being scored, so combining pulls the comparison point up exactly when,
and only when, this contig's markers genuinely disagree with the rest.

**Composition/coverage outliers** (`src/model/outliers.js`): per-bin centroid distance expressed as a z-score
(standard deviations above the bin's own mean distance), not a raw distance — makes a naturally tight bin and a
naturally diffuse one comparable on the same scale. Composition and coverage z-scores are combined via `max()`,
not average, deliberately: either axis alone being unusual is real signal, and averaging would dilute a contig
that's a strong outlier on just one axis. Coverage is only used when *every* contig in the bin has depth data
(brief's optional input); otherwise the combined score falls back to composition alone rather than partially
scoring a mix.

**Coverage table parser** finally implemented (`src/parsers/coverage-table.js`), including MetaBAT2's real
`jgi_summarize_bam_contig_depths` shape (`contigName/contigLen/totalAvgDepth/sample.bam/sample.bam-var/...` —
interleaved depth/variance columns), not just a generic numeric table. `sniff.js` disambiguates it from a
2-column contig->bin table using the column-count heuristic already documented back in Phase 1 (2 columns
defaults to contig-bin-table; 3+ columns, or a coverage-shaped header, is a coverage table).

**Kraken2 per-contig parser, a documented scope decision, not the brief's literal wording**: the brief says this
input "reuses the eDNA Explorer's Kraken2 report parser" — but that parser (`breport.js`) reads Kraken2's
*aggregate* `.breport` (one row per taxon across the whole run, no per-query field at all), which structurally
cannot answer "what did Kraken2 call this specific contig." What actually carries a per-query call is Kraken2's
other standard output format (`--output`, not `--report`): one row per query, `C/U`, seq ID, taxID, length, LCA
string. `src/parsers/kraken2-contigs.js` parses that instead, with the scope note recorded in the file itself,
not just here. The taxonomic-disagreement check built on it is deliberately simple — exact-taxID majority-vote
mismatch within a bin, not lineage-aware — since a real lineage-aware version would need a tree covering
whatever arbitrary taxa Kraken2's full database can call, unlike the marker-gene provenance check (which only
ever touches the ~5,000 taxa this app's own fixed 40-family reference set references). Flagged as a possible
future upgrade if a student also loads a full-run `.breport` alongside the per-contig file, not built now.

**Per-contig marker contribution** (`src/model/bin-summary.js`'s `computeMarkerContributions`): straightforward
given Phase 3's per-contig tags already exist — unique (found on no other contig in the bin) vs. redundant
(found elsewhere too), the exact "low-risk removal candidate" signal the brief describes.

**UI**: one new "Outlier & disagreement flagging" card, not five, sorted by a `flagCount` (how many of the five
signals are "hot" for that contig) then by combined composition/coverage z-score. Uses the reconciled putative
MAG's full contig set (core + disputed) as "the bin" for centroid/marker-contribution purposes when 2+ tools are
loaded, or the single table's own bins otherwise — the best available guess at "this genome" either way. Only
rendered when at least one bin table is loaded, since every signal here is relative to a bin.

**A real bug caught by browser testing, not the unit tests**: an empty-bin edge case (a putative MAG whose only
contig has no marker calls at all) made `computeBinTaxonomicConsistency` return an empty `Map`, so
`.get(contigId)` came back `undefined` rather than the function's documented `null` — and the render code's
`=== null` check missed it, printing the literal string "undefined" in the table. Every one of
`marker-taxonomy.test.js`'s unit tests happened to exercise bins with at least one contig lacking a call, but
never one where the *entire* group has none — the exact gap between unit coverage and an end-to-end run with
real (sparse, uneven) demo data. Fixed with `?? null` at the call site.

**Tests**: 32 new (`coverage-table.test.js`, `kraken2-contigs.test.js`, `outliers.test.js`,
`marker-taxonomy.test.js`, plus additions to `bin-summary.test.js` and `sniff.test.js`) — 109 total, all passing.
Verified in-browser with new example coverage (`examples/marker-gene-demo.coverage.tsv`) and Kraken2
(`examples/marker-gene-demo.kraken2.tsv`) files loaded alongside both existing bin tables: real NCBI-derived
taxonomic distances computed correctly (the two demo marker fixtures are genuinely unrelated organisms, so a
distance of 8 ranks is real signal, not a bug), all five signal columns populated or correctly marked n/a, no
console errors.

## Phase 7 findings — interactive reassignment

**Scope decision, stated up front rather than silently simplified**: a rectangular drag-select, not a freeform
lasso. A true lasso is arbitrary-polygon point-in-polygon hit testing — a real feature in its own right, not
just a naming difference — and a rectangle covers the same "select a visual cluster" need for a teaching tool's
scatter plot at a fraction of the implementation and testing cost. Recorded in `scatter-geometry.js`'s header,
not just here.

**Split into testable pure logic vs. DOM glue, the same pattern Phase 3's search code and Phase 5's matching
logic already used**: `src/model/working-assignment.js` (state transitions: derive the starting assignment,
reassign a set of contigs to a bin, merge, generate a fresh bin ID — all plain `Map` in, `Map` out, no DOM) and
`src/viz/scatter-geometry.js` (coordinate scaling and rectangle hit-testing, no DOM) are both fully unit tested;
`src/viz/scatter.js` is a thin, deliberately minimal layer on top that only builds SVG elements and wires mouse
events — genuinely hard to unit test meaningfully in Node without a real DOM, so it wasn't, and instead got
exercised directly in the browser (see below).

**Everything reduces to one primitive**: move a contig, split a bin, merge two bins, and spin a new bin out of a
flagged cluster all turn out to be the exact same operation — reassign a set of contig IDs to a target bin ID,
existing or brand new. `working-assignment.js` only exports that one primitive (`reassignContigs`) plus thin
convenience wrappers (`mergeBins` is just `reassignContigs` over one bin's current members), rather than four
separate code paths that would all need to stay in sync.

**Starting point for the editable assignment**: the brief doesn't say what a student edits *from* when multiple
tools were loaded. Decided: the Phase 5 reconciliation's majority vote (a contig with no majority at all — e.g.
every tool split evenly — starts unassigned rather than arbitrarily picked for it) when 2+ tools are loaded, or
the single table's own bins otherwise. No bin table loaded at all means nothing to reassign into, so the whole
section is hidden rather than shown empty — matches the existing `tools.length > 0` gate the outlier-flagging
card (Phase 6) already uses for the same reason.

**Live recalculation reuses Phase 4's `computeBinSummaries` unchanged** — the working assignment is converted
back to the same `{contigId, binId}[]` shape via `assignmentToRows` and fed straight through the existing
completeness/redundancy/N50/tier pipeline. No new "live stats" logic was needed; the existing model functions
were already pure enough to call again on demand.

**`beforeunload` guard**: a module-level `hasUnsavedReassignments` flag, set only inside the one reassignment
primitive's call site, checked (not re-registered) by a single `beforeunload` listener attached once at load.
Verified directly in the browser, not just by code inspection, since a synthetic `beforeunload` dispatch is easy
to get subtly wrong: dispatching before any reassignment confirmed `defaultPrevented` stays `false` (no nag for
a student who's only exploring, per the brief's explicit requirement), and after one programmatic drag-select +
move, confirmed `defaultPrevented` flips to `true`.

**Verified end-to-end in the browser** (not just unit tests, since the DOM-wiring layer has none): loaded the
demo assembly with both example bin tables (2 tools -> reconciliation -> interactive section appears), performed
a real rectangular drag over the rendered SVG using coordinate-mapped mouse events, confirmed the selection
count and target-bin dropdown populated correctly, moved the selection to a different MAG via the dropdown +
button, and confirmed the working-bins table recalculated immediately with the correct new contig counts and
total lengths — matching hand-computed expectations from the underlying contig lengths. Axis-switching
(re-`setDataPoints` with a different scalar) re-rendered all 3 points correctly. No console errors throughout.

## Phase 8 findings — comparison views and QC across the full set

**Scope decision on the redundancy check's similarity signal**: composition (mean tetranucleotide vector) and
mean GC, not marker-family overlap, decide whether two putative MAGs look like the same organism split across
bins. Marker-family overlap is the right signal for whole-*bin* completeness/redundancy (`bin-summary.js`) but is
useless here: the 40-family reference set is universal single-copy genes, so two MAGs from unrelated species will
both hit most of the same families at similar completeness — family identity carries no organism-identity signal
once you're comparing across bins rather than within one. Composition is the same signal `outliers.js` already
uses to separate one organism's contigs from another's inside a bin; reusing it here (in a new
`src/model/mag-redundancy.js`, since the shape of the computation — pairwise across MAGs, not per-contig within
one bin — doesn't fit `outliers.js`'s existing function signatures) keeps one consistent "what tells organisms
apart" answer across the app rather than inventing a second one.

**Completeness/contamination-across-all-MAGs scatter needed no new statistics** — it's `bin-summary.js`'s
existing `computeBinSummaries` called against a pseudo-assignment built from each putative MAG's
`coreContigIds`/`disputedContigIds` (Phase 5's `reconcileBins` output) relabelled by `magId` instead of a tool's
own bin id. `src/model/qc-comparison.js` only picks and reshapes for display (`buildMagQcPoints` for the scatter
points, `pickComparisonBins` for the good/bad choice) — no new completeness/redundancy math, mirroring the
Phase 7 lesson that `computeBinSummaries` is already generic enough to call again on any assignment shape.

**Good-vs-bad bin choice**: `completeness - redundancy` score across bins with 2+ contigs (a single-contig bin
has no cluster shape to visually contrast). Ties/insufficient data (fewer than 2 eligible bins — common with the
small example dataset, which has only 1 core contig per MAG) surface as an explicit "need at least two..." message
rather than a broken or empty plot; verified this exact fallback text renders correctly in-browser against the
shipped example data before testing the populated path.

**Verified in the browser** with two paths, since the shipped example dataset (`examples/marker-gene-demo.*`) is
too small to exercise the good/bad comparison or trigger a redundancy flag (each MAG has only 1 core contig):
(1) loaded the real example FASTA + both bin tables + coverage table via a synthetic `DataTransfer` file-input
event (no real file dialog available to automation) — confirmed the QC card renders, the completeness/redundancy
scatter plots both real MAGs at the correct SVG coordinates, and the good/bad section correctly falls back to its
"need at least two putative MAGs with 2+ contigs" message rather than erroring. (2) Called
`computeBinSummaries`/`pickComparisonBins`/`computeMagRedundancy` directly in the loaded page's console against a
synthetic 3-MAG scenario (two MAGs with near-identical composition, one clearly different) — confirmed
`pickComparisonBins` picked the right good/bad MAGs by score, and `computeMagRedundancy` flagged exactly the
near-identical pair (`compositionSimilarity: 1, likelyDuplicate: true`) while correctly not flagging either
against the distinct third MAG (`compositionSimilarity: 0.69, likelyDuplicate: false`). No console errors in
either path. All 134 tests pass (8 new, across `test/mag-redundancy.test.js` and `test/qc-comparison.test.js`).

## Phase 9 findings — export and site chrome

**About/FAQ, header/footer layout, and the responsive breakpoint already existed** in `index.html`/`styles/main.css`
before this session — built early (matching the eDNA Explorer's chrome) rather than deferred to last, unlike the
brief's phase ordering. Phase 9's actual remaining scope was the three export outputs only; no site-chrome work
was needed.

**FASTA extraction reuses the original bytes verbatim rather than reformatting.** `fasta-index.js`'s `faiEntry`
already records each contig's first-sequence-line byte offset plus its line width/wrapping — so
`computeSequenceByteSpan` (`src/model/fasta-extract.js`) only needs to compute *how many bytes* of the original
file cover a contig's sequence lines (accounting for a possibly-shorter last line), then that exact byte range —
newlines and all — is sliced straight out and reused under a new header. No re-wrapping logic was needed because
the original wrapping is already correct; the arithmetic was verified by hand for four cases (exactly one full
line, several full lines plus a short remainder, an exact multiple of the line width, and CRLF line endings)
before trusting the browser test, following the Phase 5/6 lesson above about checking metric arithmetic by hand
rather than trusting a first-draft test.

**Non-uniform-wrapping contigs are skipped from extraction, not guessed at.** `faiEntry.uniform === false` means
the file's line lengths weren't consistent for that contig, so there's no single "bytes per line" to compute a
byte span from without literally re-scanning the file — which would defeat the whole point of slicing off the
existing index. `computeSequenceByteSpan` returns `null` in that case; `planFastaExtraction` collects those as
`skippedContigIds` per group, and the UI surfaces an explicit "N contig(s) skipped (non-uniform FASTA line
wrapping...)" note rather than silently producing a truncated or corrupt sequence.

**Export always reads the live `workingAssignment`, not the originally loaded bin tables** — same pattern as
Phase 7's live-recalculating stats table, and matches the brief's explicit "revised... table, reflecting any
manual reassignment made in the session" wording. No separate "export original" option was built since nothing in
the brief asked for one.

**Verified end-to-end in the browser** against the real example dataset (`examples/marker-gene-demo.*`, loaded via
a synthetic `DataTransfer` file-input event): (1) clicked the assignment-CSV and bin-summary-CSV export buttons
with `URL.createObjectURL` monkey-patched to capture the generated Blobs instead of letting the browser save them
— confirmed the assignment CSV listed all 3 contigs against their reconciled MAGs, and the summary CSV's numbers
matched the on-page reconciliation table exactly. (2) Clicked "Download bin FASTA" for `MAG_1` — confirmed the
extracted multi-FASTA had both of `MAG_1`'s contigs under their original headers with original line wrapping
intact, zero skipped contigs, and (checked by hand against the raw source file's own sequences) the extracted
sequence lengths (1889 + 1581 = 3470 bp) matched `MAG_1`'s own reported total length exactly. (3) Resized the
viewport to a 375×812 mobile preset and reloaded the same data — confirmed the assembly summary, reconciliation
table, and the new export card all render legibly stacked single-column with wide tables scrolling inside their
own container rather than the page. No console errors in any path. All 145 tests pass (9 new, across
`test/fasta-extract.test.js` and `test/export-csv.test.js`).

This closes out all 9 phases of the brief's plan. See `docs/HANDOVER.md`'s "Open items still worth doing" for
what's left (real-data calibration follow-ups, not scoped phase work).

**Tests**: 17 new (`working-assignment.test.js`, `scatter-geometry.test.js`) — 126 total, all passing.
