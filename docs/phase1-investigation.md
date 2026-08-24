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
improved 3.17s → 2.95s (~7%) on the same 50MB benchmark — translation was no longer the dominant cost after the
first optimization round, so other steps (line-by-line parsing, `seq +=` accumulation, the composition scan, GC
counting) now make up a larger share of the total. Reported honestly rather than implying a bigger win than the
full-pipeline number shows; a further round on those other steps would be the next place to look if more
throughput is needed later.

## Still open / needs your call before I start writing code

1. **Reduced alphabet scheme** — I'll default to a published one (e.g. Murphy10) unless you have a preference.
2. **Clustering method for step 1** — greedy approximate (minhash/k-mer-profile) vs. exact — I'd default to approximate for speed, willing to revisit if it under- or over-clusters during calibration.
3. **taxdump merged/deleted taxID handling** — default plan above (map merged forward, warn-and-drop truly dead ones) — confirm this is acceptable rather than something that should hard-fail the build.
4. Ready to proceed to repo scaffolding (site skeleton + `build/` stubs) once you've reviewed this, unless you'd rather adjust something above first.
