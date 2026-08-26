# SCG marker-gene recall — BLAST-based external verification

Status: investigation complete, two threshold changes shipped (see "Changes shipped"). Two follow-up questions
about closing the remaining gap were also investigated and **not adopted** — see "Gapped extension" and
"Seed-count-only classification (no extension at all)" below.

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
protein — not another threshold nudge.

Two structural alternatives were prototyped and tested against the same BLAST ground truth as follow-ups.
Neither was adopted, but both are worth recording so the reasoning doesn't have to be redone.

## Gapped extension (prototyped, not adopted)

An indel between query and reference throws off `extendUngapped`'s register permanently — every residue after
it compares against the wrong reference position, scores as noise, and X-drop halts immediately, no matter how
well the alignment resumes a few residues later in the correct frame. `minCoverage` only controls how much of
that truncated alignment is required, not the truncation itself.

**Implemented** `extendGapped` in `src/model/marker-genes.js`: a banded, X-drop-bounded dynamic-programming
extension with a linear gap penalty (deliberately linear rather than affine open+extend, for a first prototype),
same left/right-from-anchor structure as `extendUngapped`. Wired in behind a `useGappedExtension` param,
**off by default** — shipped behavior is unaffected unless a caller opts in. Sanity-checked two ways: with gaps
effectively disabled (a very high gap penalty) it reproduces `extendUngapped`'s score exactly; on a synthetic
3-residue deletion, ungapped extension covers 31% of the reference before stopping, gapped extension crosses it
and covers 100%, tripling the score.

**Re-verified against the same BLAST ground truth**, stacked on top of the shipped `minCoverage: 0.3`:

| | Contigs called | Families | Recall (exact) | BLAST-unsupported | Search time (full assembly) |
|---|---|---|---|---|---|
| Ungapped (shipped) | 1,133 | 36 | 74.2% | 6.7% | 8.4s |
| Gapped, bandWidth=16 | 1,172 | 37 ⚠ | 77.1% | 7.6% | 36.5s |
| Gapped, bandWidth=4 (narrowest tested) | 1,165 | 37 ⚠ | 76.4% | — | 17.9s |

A real, independent ~3-point recall gain — but two problems that ruled it out for now:

- **The runtime floor doesn't move much with band width.** Even at the narrowest band tested (4, vs. the
  default 16), search time was still ~2.1x the ungapped baseline for essentially the same recall. The overhead
  is inherent to running a DP table instead of a linear walk, not something a band-width knob can tune away —
  the realistic choice is "roughly 2-4x slower" vs. "don't," not a dial with a cheap end.
- **A new spurious family appears at every band width tested** (COG0096 — not supported by BLAST at all), the
  same failure shape as pushing `minCoverage` too low. Fixing it needs its own calibration pass (the gap penalty
  is the suspect, not band width, since it persists even at the narrowest band) — on top of the fact that every
  existing threshold (`minScore`/`minCoverage`/`minMargin`/`minRepresentatives`) was tuned against ungapped
  score distributions and would need re-deriving for gapped ones before this could ship.

Given this tool's explicit "lightweight browser visualizer, not production" framing (the same framing behind
`minSeedStride`/`maxSeedsPerSegment`'s careful seed-volume budgeting elsewhere in `marker-genes.js`), a
consistent 2-4x slowdown of the search stage for +3 points of recall, bundled with a new precision problem that
would need its own calibration investment, wasn't judged worth pursuing further. The code stays in the module
as a tested, off-by-default option — reachable via `{ useGappedExtension: true }` — in case that judgment
changes later (a bigger runtime budget, or recall becoming a higher priority).

## Seed-count-only classification (no extension at all)

The opposite direction: what if extension were dropped entirely, and a family call was based purely on how many
independent seeds land on a diagonal, with a higher count cutoff instead of a score/coverage bar? Extension is
by far the most expensive per-candidate step, so this was worth checking even though the expectation going in
was that it would trade precision away.

Swept the raw seed-count cutoff (using the same seeding/bloom/diagonal-counting path production uses, just
stopping before any `extendUngapped`/`extendGapped` call) against the same BLAST ground truth:

| min seed count | Contigs called | Recall (exact) | BLAST-unsupported |
|---|---|---|---|
| 2 | 21,852 | 79.9% | 94.2% |
| 4 | 12,670 | 79.4% | 90.0% |
| 6 | 1,721 | 76.6% | 33.0% |
| 10 | 1,274 | 74.5% | 11.9% |
| 20 (highest tested) | 1,191 | 72.7% | 8.1% |
| *shipped (extension-based)* | *1,133* | *74.2%* | *6.7%* |

The expectation held, decisively. At any cutoff loose enough to be useful (≤4 seeds), 90%+ of calls are false
positives — a raw seed count can't distinguish "several short exact matches that are part of one real
homologous alignment" from "this particular short reduced-alphabet motif just happens to recur across many
unrelated diagonals by chance." Extension is exactly the step that makes that distinction, since a coincidence
can share a few scattered k-mers with something but can't fake a long run of good BLOSUM62-scored residue
similarity between them. Even at the highest cutoff tested (20), the result is worse than the shipped approach
on **both** recall (72.7% vs. 74.2%) and precision (8.1% vs. 6.7% unsupported) — there's no cutoff in this sweep
that matches shipped recall while beating shipped precision. Extension isn't overhead that could be trimmed
away for a faster classifier; it's the step doing the actual discrimination. Not pursued further — no code
changed as a result of this check (a throwaway script, not added to the module).

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
