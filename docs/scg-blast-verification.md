# SCG marker-gene recall — BLAST-based external verification

Status: investigation complete, two threshold changes shipped as a result (see "Changes shipped" below).

## Why this is a different check from Phase 3's calibration

[phase1-investigation.md](phase1-investigation.md)'s "Phase 3 calibration findings" measured the search
against **held-out sequences from the same reference collection** the shipped index was built from — a
trustworthy *relative* signal between families, but explicitly flagged there as not a substitute for the
brief's own ask: testing against real assemblies not specifically enriched for these 40 families. This
investigation is that missing piece — an independent, real-world recall check using an actual BLAST-family tool
(DIAMOND) as ground truth, run against a real 45,388-contig, ~195Mb multi-tool-binned assembly (5 binning
tools: MaxBin2, MetaBAT2, SemiBin2, VAMB, CONCOCT) supplied for testing, not part of the reference set at all.

## Method

1. Built a DIAMOND protein database from **`reference-data/scg40_raw.fasta`** directly — the full 69,208-sequence
   raw reference set (pre-clustering), not the shipped 39,854-sequence clustered index, so the BLAST side isn't
   subject to the same clustering choices as the k-mer side.
2. Ran `diamond blastx` (v0.8.38; BLAST+ 2.15.0 is equivalent but much slower at this scale) of the real
   assembly's 45,388 contigs against that database: `--evalue 1e-10 --max-target-seqs 5`.
3. Filtered hits to `pident >= 30` and alignment `length >= 100` aa — a permissive but standard protein-homology
   confidence bar — then took each contig's single best hit by bitscore. This produced a **ground-truth set of
   1,386 contigs BLAST calls as carrying a marker gene**, spanning 36 of the 40 families (the other 4 are
   presumably genuinely absent or too divergent from this reference set for BLAST to catch either, at this
   e-value/identity bar).
4. Ran the actual production search (`src/model/marker-genes.js`'s `searchContigForMarkers`, the same code path
   the app uses) over the same assembly, and compared per-contig family calls against the BLAST set.

Family-level result, before any changes: **both approaches find exactly the same 36 of 40 families** — the
built-in search isn't missing any family that's genuinely present, only under-calling *which contigs* carry it.

## Baseline: 61.9% contig-level recall

Of the 1,386 BLAST-positive contigs, the k-mer search (default thresholds) matched the exact family on 858
(61.9%) and missed 491 entirely. Precision was high — of the search's 948 total calls, only 53 (5.6%) had no
BLAST support at all — confirming the module's own documented framing ("under-calls divergent lineages before
it over-calls").

## Root-cause tracing of the misses

Walked the 491 misses through each stage of the search pipeline in order, using
`computeFamilyCandidates`/`extendUngapped`/`lookupPopulatedIndex` directly (not just the final yes/no), to find
where each one actually dropped out:

1. **Sparse-key pruning** (`MIN_HITS_PER_KEY=3` in `build/02-index.js` — any reduced-alphabet 9-mer occurring
   fewer than 3 times across the reference set is dropped from the shipped index entirely): rebuilt the raw,
   unfiltered key-occurrence map over the exact clustered reference set (`build/01-cluster.js`'s
   `scg40_clustered.fasta`, 39,854 sequences) and checked whether each miss's only shared k-mer with the true
   family had been pruned. **Only 2 of 451** "no candidate at all" misses were explained by this.
2. **Hit-list capping** (`MAX_HITS_PER_KEY=20` — a key that's kept can still only store an evenly-strided sample
   of its hits): checked the actual shipped index's stored hit lists directly. Same **2 of 451**.
3. **Two-hits-per-diagonal gate** (`MIN_SEEDS_PER_DIAGONAL=2` — extension is never attempted from a single,
   unconfirmed seed): **27 of 451** had fewer than 2 independent seeds landing on the correct diagonal against
   the true reference family.
4. **The actual dominant cause — coverage of an ungapped extension**: the remaining **424 of 451 (94%)** *did*
   seed and *did* extend, often to a strong BLOSUM62 score (some over 1000, against `minScore=50`), but failed
   `minCoverage=0.5` — coverage is measured as aligned length over the *full reference protein length*, and
   extension is strictly ungapped. A genuine 30–70%-identity homolog's indels stop an ungapped X-drop extension
   well short of 50% of a 200–1500-residue reference protein, regardless of how good the local match is. This is
   an architectural interaction (ungapped extension + full-length coverage denominator), not a bug or a sign of
   a bad index.

So the initial hypothesis going in — that rare k-mers pruned during index-building were the main cost — turned
out to explain under 1% of the gap. The real lever was `minCoverage`.

## Changes shipped

Both landed in `data/scg40-thresholds.json` (the same file `build/03-calibrate.js` populates; these two
overrides came from this external-verification investigation, not from that script — see "Relationship to
Phase 3 calibration" below).

**1. `minRepresentatives: 1` override for 9 families** (COG0012, COG0016, COG0018, COG0085, COG0124, COG0201,
COG0202, COG0215, COG0552), alongside the two Phase 3 already found (COG0495, COG0525). 38 of the 491 misses
had the correct family found with a strong, unambiguous score (several >1000 bitscore, >95% identity) but
rejected purely because only 1 reference representative cleared threshold against the default of 2. Recovers
exactly those 38 contigs, no measurable precision cost (53 unsupported calls before and after).

**2. Global `minCoverage: 0.5 → 0.3`**. Tested the full sweep before choosing this value:

| minCoverage | Contigs called | Families found | Recall (exact family) | BLAST-unsupported |
|---|---|---|---|---|
| 0.5 (previous default) | 986 | 36 | 64.6% | 5.4% |
| 0.4 | 1,061 | 36 | 69.6% | 5.9% |
| **0.3 (shipped)** | **1,133** | **36** | **74.2%** | **6.7%** |
| 0.2 | 1,217 | 37 ⚠ | 78.4% | 8.5% |
| 0.15 | 1,266 | 38 ⚠ | 80.2% | 10.2% |
| 0.1 | 1,361 | 38 ⚠ | 81.9% | 14.6% |

0.3 is the last value that still finds exactly the same 36 families BLAST does — 0.2 and below start calling
families (e.g. COG0091 at 0.2) that BLAST doesn't support at all, a genuine spurious family-level call, not just
a borderline contig — and it sits at the elbow of the recall/precision curve (+9.6 points of recall for +1.3
points of unsupported calls; returns diminish sharply below it).

**Combined effect**: contig-level recall (exact family) moved from **61.9% → 74.2%**, BLAST-unsupported calls
from 5.6% → 6.7%, family-level detection unchanged at 36/40.

## What's left (a real ceiling, not a threshold to tune further)

After both changes, the remaining gap is almost entirely the ungapped-extension/coverage interaction described
above. Closing more of it would need a structural change — gapped extension, or a coverage threshold that
scales with how much of the contig's own length is available rather than a flat fraction of the full reference
protein — not another threshold nudge. Worth revisiting if recall becomes a priority again, but out of scope
for this pass.

## Relationship to Phase 3 calibration

`build/03-calibrate.js` derives its overrides from held-out members of the *same* reference collection the
index came from — a relative, internal signal (see phase1-investigation.md's caveat on this). This
investigation's overrides came from an independent, real-assembly, real-BLAST source instead. Re-running
`build/03-calibrate.js` in the future (e.g. after a reference-set update) will **not** reproduce or preserve the
9 `minRepresentatives` overrides or the `minCoverage` change added here — it only ever tightens/loosens what its
own held-out data supports, and doesn't currently look at `minCoverage` at all (flagged as out of scope in its
own header comment, for the same reference-internal-diversity overfitting reason). Anyone re-running that script
should carry these two changes forward by hand, or fold this verification's method into it, rather than assume
its output supersedes this file.

## Reproducing this

Not checked into the repo as a script (this was an ad-hoc investigation, not a maintained tool) — the recipe,
if repeated against a new real assembly:

1. `diamond makedb --in reference-data/scg40_raw.fasta -d scg40_db`
2. `diamond blastx -q <assembly.fasta> -d scg40_db.dmnd -o hits.tsv -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp --evalue 1e-10 --max-target-seqs 5`
3. Filter to `pident >= 30 && length >= 100`, take best hit per contig by bitscore, family = `sseqid` up to the
   first `.`.
4. Run `searchContigForMarkers` (via `src/model/marker-genes.js`, translating with `src/model/translate.js`,
   same pattern as `test/bench-marker-genes.js`) over the same assembly and compare per-contig family calls.
5. For root-causing specific misses, `computeFamilyCandidates` (exported from `marker-genes.js`) exposes the
   pre-margin/pre-agreement candidate list — the same entry point `build/03-calibrate.js` uses — and
   `lookupPopulatedIndex`/`extendUngapped` are also exported for stepping through the seed/extend pipeline by hand.
