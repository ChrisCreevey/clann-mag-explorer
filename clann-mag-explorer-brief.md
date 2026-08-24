# Clann MAG Explorer — design brief for Claude Code

## Purpose

This document specifies a browser-only web app for exploring and refining prokaryotic metagenome-assembled genome (MAG) binning results from a single assembly. It is the fifth tool in the Clann suite, alongside [Clann Tree Viewer](https://chriscreevey.github.io/clann-tree-viewer/), [Clann BLAST Explorer](https://chriscreevey.github.io/clann-blast-explorer/), [Clann Pangenome Explorer](https://chriscreevey.github.io/clann-pangenome-explorer/), and the eDNA Explorer. It matches those in architecture and aesthetic: vanilla HTML/CSS/JavaScript, no build step, no external dependencies, GPL-2.0 licence, theme-aware CSS.

The tool does not run any binning algorithm itself. It takes an assembly and one or more independent contig→bin assignments (from different binning tools run on that same assembly), computes what it can from the sequences directly, and gives students an interactive way to see where those binning tools agree, where they disagree, and to manually reassign contigs in response. It sits closer to Anvi'o's manual refinement interface than to a QC dashboard, but working from finished binning predictions rather than from an assembly with no bins at all.

Primary audience: undergraduate and postgraduate students learning genome-resolved metagenomics. Secondary audience: researchers doing a first-pass reconciliation of multiple binning tools before committing to a formal consensus method such as DAS_Tool.

## Positioning within the Clann suite

Working name: **Clann MAG Explorer**. Suggested repository: `clann-mag-explorer`, hosted at `chriscreevey.github.io/clann-mag-explorer/`.

Cross-links belong in the header/footer and About section, exactly as the other tools link to each other. A specific contig or gene worth checking can cross-link to [Clann BLAST Explorer](https://chriscreevey.github.io/clann-blast-explorer/); a finished set of high-quality MAGs is a plausible input to [Clann Pangenome Explorer](https://chriscreevey.github.io/clann-pangenome-explorer/) for comparative analysis.

**Prerequisite before build starts**: the marker-gene identification module (see below) depends on a reference set Chris will add directly to the repository — this is not something Claude Code needs to source itself. See "Provisioning the marker-gene reference set" for what's expected and how it's used.

## Background: the problem this tool addresses

Automated binning tools (MetaBAT2, CONCOCT, MaxBin2, and others) disagree with each other, sometimes substantially, on how to group the contigs of a single assembly into genome bins. DAS_Tool exists to reconcile this, but it works at whole-bin granularity: it scores each tool's candidate bins using single-copy marker gene counts, then greedily picks the best-scoring whole bin wherever tools' candidates overlap, discarding the rest. A winning bin can still contain contigs that other tools placed elsewhere, since DAS_Tool never inspects agreement contig by contig.

This tool works at the other granularity: for each contig, how many of the loaded binning tools agree on where it belongs. That agreement fraction, computed by matching each tool's bins to one another via contig overlap, gives a **high-confidence core** per putative genome (contigs every tool agrees on) and a **disputed set** (contigs where tools disagree), without needing marker-gene search. Combined with statistical outlier detection within a bin (a contig sitting far from its bin's composition/coverage centroid) and, optionally, a per-contig taxonomic call, this produces a ranked list of contigs most worth a student's manual attention, and an interactive interface to act on that list directly.

## Inputs

| Input | Required for | Format |
|---|---|---|
| Assembly (contig FASTA) | Everything; source of all per-contig statistics | `.fasta`/`.fa`/`.fna`, optionally `.gz` |
| Contig→bin assignment table (one or more) | Bin summaries; required in at least one copy, central to cross-tool reconciliation once two or more are loaded | Two-column tab/CSV: contig ID, bin ID, one file per binning tool |
| Coverage depth table | Coverage-aware outlier detection and bin coverage profiles | Tab/CSV: contig ID, per-sample depth columns |
| Per-contig taxonomic call | Taxonomic-disagreement flagging within a bin | Reuses the eDNA Explorer's Kraken2 report parser, run against the assembly's contigs rather than reads |
| Pre-computed per-contig marker-gene hits | Alternative to the built-in marker-gene module, preferred over it when present, since real profile-HMM output is more trustworthy than the browser's approximation | Investigation needed (see Phase 1) into a realistic minimal format a lightweight HMMER run could produce |

The marker-gene reference set is **not** a student-facing input — it ships as a bundled asset with the site itself. See "Provisioning the marker-gene reference set" below.

All file inputs also accept `.gz` compression, decompressed in the browser. Detection is content-based, not by file name or extension, consistent with the eDNA Explorer.

## Provisioning the marker-gene reference set

**This is a build-time asset the repository owner supplies, not something the deployed tool ever asks a student for.** Chris will add the raw reference multi-FASTA (40 gene families, amino acid sequences, header format `>COGID.fa.<id>.<locus_tag>`, ~69,000 sequences across all families) to the repository at a path agreed during setup (e.g. `/reference-data/scg40_raw.fasta` — exact path is Claude Code's call). It should be treated as present in the repo from the start of the build, not as something to prompt for later.

That raw file is the **input to Phase 1's offline build step**, not something the deployed site ships or loads at runtime, and neither is a plain reduced FASTA. That offline step does three things in sequence: pre-clustering (collapsing near-identical representatives within each family, see the marker-gene identification module below for why), building the seed index itself from the clustered result (the k-mer-to-(reference sequence, position) lookup structure the matching pipeline queries), and resolving each distinct taxID present in the reference headers against a full NCBI taxonomy dump (downloaded at build time only) to produce a compact taxID→lineage table for the marker-gene taxonomic consistency check (see below). All three are pure, deterministic computation over fixed reference data with no session-specific input, so all three belong at build time, not in the student's browser. What actually ships to the browser is the finished index, the reference sequences needed for extension scoring, and the taxID→lineage table, all pre-built, in formats the browser loads directly rather than constructs. The raw file and the build script belong in the repository so these assets can be regenerated if the clustering threshold, k-mer parameters, or taxonomy version are revisited later, but nothing downstream of that script runs client-side.

Practical consequence for the build: this needs a Node (or similar) build-time script, run outside the browser, separate from the zero-dependency constraint that governs the deployed site itself. That constraint is about what ships to the student's browser, not about how the repository's own build tooling works.

## Marker-gene identification module

This is a genuine scope expansion beyond "visualise and assess what's supplied," and is treated as such: a clearly separated, optional module rather than folded silently into the core phases. Without it, completeness/contamination simply stay "not available" unless a student supplies pre-computed marker-gene hits from a real tool. With it, every student gets an estimate with no extra pipeline step, and — because the search runs once, in the same streaming pass as everything else — bin-level completeness/redundancy recalculates live during interactive reassignment, which a pre-computed report can never do.

**The search is a single pass, not an ongoing cost.** Each contig is translated in all six reading frames (three forward, three reverse-complement) as one continuous amino-acid string per frame, roughly contig length ÷ 3, with stop codons written as a sentinel character rather than used to pre-cut the string into segments. This is searched against the marker-gene reference set once, during initial loading. The output is a tag per contig: which of the 40 gene families it hit, and how many times. Everything downstream — bin-level completeness, redundancy, and the per-contig contribution view below — is aggregation over these tags, cheap enough to redo on every reassignment.

No separate ORF-calling or segmentation step is needed. Framing each frame as one continuous string with an unmatchable sentinel at stop-codon positions makes segmentation happen implicitly, as a side effect of normal seeding and extension, rather than as up-front work: no k-mer lookup key can span a sentinel, so seeds only ever form within real coding stretches; extension halts on its own the moment it reaches a sentinel or the string's end, whichever comes first, with no separate boundary logic and no special case needed for a marker gene fragment truncated by the contig's edge rather than by a stop codon. There's also no start-codon requirement, since a fragment truncated at a contig edge won't have one in view and would otherwise be wrongly discarded.

**Matching pipeline:**
1. **Reduced-alphabet index, precomputed offline and shipped as a static asset, not built by the browser.** During the build step described above, every reference sequence is translated into a reduced amino-acid alphabet (roughly 10–12 groups of biochemically similar residues, collapsing most conservative substitutions to the same symbol), broken into overlapping fixed-length windows (k=5–6 residues), and each window's contents used purely as a hash-table lookup key, mapped to every (reference sequence ID, position) where it occurs. This is k-mer-as-lookup-key, a search-acceleration device, not a k-mer frequency profile: nothing is counted, no similarity is measured at this stage, it only answers "which reference sequences share this exact short reduced-alphabet stretch, and where." (The tool's contig-binning signature elsewhere uses actual k-mer frequency profiles for a genuinely different, compositional purpose; this step deliberately isn't that, since composition can't distinguish true homology from coincidentally similar amino acid usage, which matters given the paralog risk below.) The finished index — not the reference FASTA, not the reduced-alphabet translation step — is what the student's browser downloads and loads at page start, in a compact binary layout (parallel typed arrays, integer-encoded k-mer keys) rather than JSON, consistent with the low-memory-footprint precedent set by the Pangenome Explorer. The real (non-reduced) reference sequences ship alongside it, since extension (step 3) needs the actual residues, not just the index.
2. **Seeding per frame.** The same reduced-alphabet windowing is applied to each translated frame string, each window looked up against the shipped index. A hit returns (frame position, reference sequence ID, reference position), for every reference sequence sharing that stretch, across every family, not narrowed down at this stage.
3. **Diagonal binning and ungapped extension, scored on real amino acids.** Hits against the same reference sequence sharing a diagonal (frame position − reference position) mark a plausible ungapped alignment region. From a seed on a promising diagonal, extension proceeds left and right using the actual BLOSUM62 matrix on the real (non-reduced) residues, drawn from the shipped reference sequence data, with an X-drop stopping rule: keep extending while the running score stays within X of its best value seen so far, stop once it falls further behind than that. The reduced alphabet is used only for the seeding lookup; scoring always happens on real residues. This step runs per (frame, reference sequence) pair independently — nothing is merged or averaged across representatives at this stage, and nothing is discarded down to a single best match yet.

**Collapsing many representatives to a family-level result happens only after extension, and needs two different numbers kept apart, not one.** For a genuine marker gene, seeding typically returns hits against a large fraction of that family's ~1,700 representatives, each independently extended and scored. From that full per-representative result set:
- The **best individual score** within the family feeds the threshold and margin-against-second-best-family checks below.
- **How many representatives independently clear the threshold** feeds the multi-representative agreement check below.

Consolidating early into one profile per family (averaging or merging before scoring) would destroy the second number, since it depends on counting how many distinct representatives agree, not on a blended composite. Keeping only the single best match and discarding the rest would do the same. Both numbers are read off the same kept-apart per-representative result set, only after extension is done.

This stays computationally reasonable because the two failure directions are self-limiting in opposite, convenient ways: a non-marker stretch of translated sequence produces few or no seed hits at all, so little or no extension work follows it; a genuine marker gene does trigger extension against many representatives, but that's exactly the case where the cost is earning something, since it's the evidence the agreement check needs. An early-stop optimisation — halting extension for a family once a small number of representatives (likely 2–3, matching whatever the agreement threshold turns out to require) have independently cleared the score/coverage threshold, rather than exhaustively extending against all ~1,700 — is worth Phase 1 tuning rather than fixing now, since it trades a little confidence (seeing every representative's score) against real speed, and the right cutoff is an empirical question.

**Paralog safety, since no negative/decoy sequences exist in the reference set to calibrate against directly.** Three checks applied together against the collapsed per-family results above, not individually, given the explicit risk of mis-assigning paralogs with only positive examples to match against:
- A minimum score/coverage threshold in BLOSUM62 bit-score terms, requiring the aligned region to cover most of the shorter sequence, not just a short motif shared across unrelated families.
- **A margin requirement between the best and second-best gene family's score** — reciprocal-best-hit logic — so a hit region is only assigned where one family's best score clearly beats every other family's, not just edges it out. This is the primary defence against the marker set's own paralogous families (several of these COGs, particularly aminoacyl-tRNA synthetase-adjacent ones, have close relatives within the set).
- **Multi-representative agreement**: require more than one independent representative within the winning family to have independently cleared the threshold, not just one, made practical by the reference file's depth (~1,700 representatives per family) and the per-representative extension approach above. A genuine match clears this easily; a spurious single-sequence hit mostly won't.

A hit region failing any of the three stays unassigned. Under-calling is the accepted failure direction, consistent with the sensitivity caveat below.

**Per-contig marker contribution, the distinctive output beyond a bin-level number.** For a contig currently in a bin, distinguish:
- **Unique contribution**: marker genes on this contig found on no other contig currently in the same bin. Removing this contig loses completeness.
- **Redundant contribution**: marker genes on this contig also found on other contigs in the same bin. Removing this contig reduces redundancy (apparent contamination) without costing completeness.

A contig with high redundant and zero unique contribution is a low-risk removal candidate, made explicit and sortable rather than left for a student to infer from clustering alone. This becomes a third signal in the outlier-flagging view, alongside composition/coverage distance and cross-tool disagreement.

**Marker-gene taxonomic consistency, using provenance the search already produces for free.** Each reference sequence's header carries an NCBI taxID (`>COGID.fa.<taxID>.<locus_tag>`), so every called marker gene has not just a family assignment but a **provenance**: the taxID of its top-scoring representative, treated as "nearest known relative in the reference set," not a species-level identification. Collecting these taxIDs across all of a bin's called markers and computing their lowest common ancestor gives a chimerism signal analogous to GUNC's lineage-homogeneity check, computed here as a side effect of the marker search rather than a separate whole-genome classification pass. A bin whose markers only agree at a coarse rank (phylum or class, say, where the marker set would normally resolve much finer) is a chimerism candidate; per-contig, whichever contig's markers sit furthest from the bin's consensus lineage becomes a fourth outlier signal, alongside composition/coverage distance, cross-tool disagreement, and marker redundancy.

This needs a taxID→lineage lookup that doesn't exist in the reference FASTA itself, unlike the eDNA Explorer, where the input report already encodes the full taxonomy tree via indentation. That lookup is built the same way as the marker-gene index: offline, once, during the build step (see "Provisioning the marker-gene reference set" above) — the build script resolves each of the reference set's few thousand distinct taxIDs against a full NCBI taxonomy dump (downloaded at build time, never shipped) and produces a compact taxID→lineage table covering only those taxa, shipped as a further static asset alongside the index and reference sequences. The tree-traversal/LCA logic itself should be **copied from the eDNA Explorer's codebase into this repo**, not imported at runtime, since these are independent static sites with no shared runtime and the no-dependency constraint rules out reaching across repos.

**Stated limitation, surfaced in the UI next to any number this module produces, not just in documentation**: this is an ungapped seed-and-extend heuristic (k-mers used only as lookup keys for seeding, never as a similarity measure) against a kept-apart, per-representative scoring scheme, not a profile HMM search. It will under-call divergent or poorly-represented lineages before it over-calls, but it is not equivalent to CheckM2 or HMMER output, and should never be presented as though it were.

**Build-time considerations, not runtime ones:**
- **Pre-clustering the reference set offline** (e.g. collapsing near-identical representatives within a family at ~90% identity) before shipping it, via a one-time build script rather than in-browser, keeps the shipped index smaller without losing the diversity the margin and consistency checks depend on.
- **Building the seed index itself offline**, immediately after pre-clustering, rather than at runtime in the browser: the index is a pure function of the reference set, with no session-specific input, so constructing it fresh in every student's browser is unnecessary repeated work. The build script's output — the finished index plus the real reference sequences needed for extension — is what ships; the browser only ever loads it.
- **Choice of shipped index format** — a compact binary layout of typed arrays rather than JSON — affects both download size and page-load parse time, and should be settled during Phase 1 rather than left implicit.
- **Empirical threshold calibration** — testing the chosen thresholds (score/coverage, family-margin, minimum agreeing-representative count, and the early-stop cutoff above) against frame translations of genomes not expected to be enriched for these 40 families — should happen offline during Phase 1, so the shipped defaults are checked rather than guessed.

## The FASTA is never held in memory

This is a hard constraint, not an implementation detail, given assemblies can run to tens of thousands of contigs and hundreds of megabytes.

- The FASTA is read once via a streaming, chunked read (`Blob`/`ReadableStream` access, not loading the whole file into a string), one contig at a time.
- On this single streaming pass, the tool computes and retains, per contig: length, GC content, GC skew, a composition signature (k-mer/tetranucleotide frequency), a coding-density estimate (six-frame translation, stop-to-stop segment lengths), the marker-gene module's tags (reusing the same six-frame translation, searched as continuous sentinel-delimited strings rather than pre-cut segments — see the marker-gene identification module below), and **a FASTA index entry**: byte offset where the sequence starts, sequence length, and line-layout (bases per line, bytes per line), the same information `samtools faidx` produces.
- The sequence itself is discarded once its per-contig record is computed. Only the numeric record and index entry are kept.
- The original `File`/`Blob` reference is kept alive for the lifetime of the session (it supports random-access `.slice(start, end)` reads without being "consumed" by the first pass), so later sequence extraction never needs the user to re-supply the file, and never needs a second streaming pass over the whole assembly. It reads exactly the bytes needed via the index, for exactly the contigs needed, as many times and in as many different MAG groupings as the session requires.
- **Known limitation, to be stated plainly rather than engineered around unless it proves worth it**: a page reload or tab close loses the live `File` reference, since browsers have no way to silently reopen a file. Computed statistics could in principle be preserved (e.g. via in-browser storage), but sequence extraction after a reload would need the user to re-supply the original file. Phase 1 should decide whether this is an acceptable limitation for a teaching tool or worth the added complexity of persistable File System Access API handles.
- Accepting a pre-existing external `.fai` index as an alternative to building one is left as an open point rather than built by default, since building one is a free byproduct of the pass the tool makes anyway, and accepting one is an extra format to validate against the FASTA it's meant to describe.

## Layout and site chrome

Matches the eDNA Explorer's layout exactly, since a student moving between the two tools should not have to relearn the interface: filters and search in the **left pane**, visualisations and tables in the **right pane**, header and footer mirroring the eDNA Explorer's structure (site title/logo linking home, navigation to the other Clann tools, GitHub link, licence and version in the footer, consistent theme-aware CSS).

**Unsaved-work protection**: since manual contig reassignment can represent a substantial amount of a student's work, and since a lost `File` reference after a reload also means lost sequence-extraction capability (see above), the tool must intercept navigation away from the page (`beforeunload`) once any reassignment has been made in the session, and show the browser's native "are you sure you want to leave" confirmation. This should arm as soon as the first reassignment happens, not before, so it doesn't nag a student who's only exploring.

## Left pane: filters and search

- Filter contigs by length, GC%, coding density
- Filter by agreement fraction (e.g. show only contigs where tools disagree)
- Filter by bin (from any loaded tool) or by "unbinned"
- Contig ID / bin ID search, highlighting matches across all open views
- Toggle which loaded binning tool(s) are active in the comparison views, so a two-tool or three-tool comparison can be explored without reloading

## Right pane: analyses and visualisations

### Per-contig properties (computed on load)
- Length, GC content, GC skew, composition signature, coding-density estimate, coverage per sample if supplied, taxonomic call if supplied

### Per-bin / per-putative-MAG summary
- Contig count, total length, N50/L50, mean GC, coverage profile across samples
- SCG-based completeness/redundancy, from the built-in marker-gene module or from supplied pre-computed hits (the latter preferred when both are present)
- Marker-gene consensus lineage (lowest common ancestor across the bin's called markers' provenance) and the rank at which it resolves, flagging bins where markers only agree at an unexpectedly coarse rank as chimerism candidates
- A MIMAG-style quality tier computed from whichever completeness/contamination figures are available, with thresholds shown and adjustable, consistent with the eDNA Explorer's transparent-calculation philosophy

### Cross-tool reconciliation (core feature)
- Bin-matching across loaded tools by contig overlap, so equivalent bins from different tools are recognised as the same putative genome despite different bin IDs
- Per-contig agreement fraction across loaded tools
- High-confidence core set per putative MAG
- Disputed contig list, ranked by how split the disagreement is
- Side-by-side view of what each tool's version of a given putative genome includes

### Outlier and disagreement flagging
- Contigs statistically distant from their bin's composition/coverage centroid
- Overlaid with the cross-tool agreement signal, so a contig that is both a statistical outlier and cross-tool disputed is flagged as the strongest reassignment candidate
- Taxonomic-disagreement flag where a contig's Kraken2 call differs from the rest of its bin
- Marker-gene contribution per contig (unique vs. redundant, see the marker-gene module above), sortable, and combined with the other signals so a contig that is an outlier, disputed, and purely redundant in its marker contribution is flagged as a very strong removal candidate
- Marker-gene taxonomic consistency per contig: how far its markers' provenance (nearest-reference taxID, see the marker-gene module above) sits from the bin's consensus lineage, a distinct signal from the Kraken2 flag above since it comes from marker-gene homology rather than whole-contig classification, and is available even when no Kraken2 run against the assembly is supplied

### Interactive reassignment
- Scatter/parallel-coordinates view (composition, coverage, GC, length) with lasso-style selection
- Move a contig between bins, split a bin, merge two, or spin a new bin out of a cluster of flagged contigs
- Every bin-level summary statistic recalculates live as reassignment happens

### Comparison and QC across the full set
- Good-bin vs. bad-bin side-by-side view (tight, well-separated cluster vs. scattered/fragmented one), for direct teaching illustration
- Completeness vs. contamination scatter across all putative MAGs, coloured by which tool(s) support each one
- Redundancy check across finished MAGs to flag likely duplicate genomes split across separate bins

## Export

- Revised contig→bin assignment table, reflecting any manual reassignment made in the session
- Per-MAG FASTA extraction, built via `Blob.slice` against the first-pass index rather than a second full file read, supporting extraction of the same contig into multiple different MAG groupings (one tool's version, another tool's version, the reconciled core set) without re-parsing
- Summary table per MAG: length, GC, N50, tool-agreement fraction, completeness tier

## Suggested build phases

1. **Investigation and data model.** Confirm streaming FASTA parsing and index-building strategy (byte offsets, line-layout, `Blob.slice` extraction) against real assemblies of realistic size. Decide the reload/File-reference limitation question. Investigate what a realistic pre-computed per-contig marker-gene hit format would look like from a lightweight HMMER run, and whether accepting an external `.fai` index is worth supporting. Design content-based detection for contig→bin tables and coverage tables. **Confirm the raw marker-gene reference file has been added to the repository (see "Provisioning the marker-gene reference set")**, then build the offline pipeline (pre-clustering, seed-index construction, and taxID→lineage resolution against an NCBI taxonomy dump, in that order), decide the shipped index's binary format, and empirically calibrate the matching thresholds (score/coverage, margin, multi-representative agreement, early-stop cutoff) against genomes not expected to be enriched for the 40 families. Port the eDNA Explorer's tree-traversal/LCA logic into this repo.
2. **Single-assembly parsing and per-contig statistics.** Streaming FASTA read, index built on the same pass, length/GC/GC-skew/composition/coding-density computed and retained, sequence discarded.
3. **Marker-gene identification module.** Load the precomputed seed index, reference sequences, and taxID→lineage table (built offline in Phase 1, shipped as static assets), run the seed-and-extend search against each contig's six-frame translation, apply the three paralog-safety checks, produce per-contig gene-family tags and provenance taxIDs during the same initial pass as Phase 2's statistics.
4. **Single binning-result loading and bin summaries.** Load one contig→bin table, per-bin summary stats (including completeness/redundancy from Phase 3's tags or from supplied pre-computed hits), MIMAG-tier calculation.
5. **Cross-tool bin matching and reconciliation.** Load multiple contig→bin tables, match equivalent bins across tools by contig overlap, compute per-contig agreement fraction, build the core/disputed contig sets.
6. **Outlier flagging and taxonomic overlay.** Composition/coverage centroid distance per bin, optional Kraken2 taxonomic-disagreement flag (reusing the eDNA Explorer's parser), per-contig marker contribution (unique vs. redundant), marker-gene provenance LCA/consensus lineage per bin and per-contig distance from it.
7. **Interactive reassignment.** Scatter/parallel-coordinates view, lasso selection, live-recalculating bin stats including completeness/redundancy, the `beforeunload` unsaved-work guard armed on first reassignment.
8. **Comparison views and QC across the set.** Good-vs-bad side-by-side, completeness/contamination scatter, redundancy check.
9. **Export and site chrome.** Revised assignment table, index-based per-MAG FASTA extraction, summary table export; header/footer/left-right layout matching the eDNA Explorer, About/FAQ, responsive layout pass.

## Open points for discussion

- Whether to support persistable File System Access API handles to survive a page reload, or accept the simpler limitation (re-supply the file after a reload) as adequate for a teaching tool.
- Whether to accept an external `.fai` index as an alternative to always building one internally.
- What a realistic, minimal pre-computed marker-gene hit input format looks like, and whether it's reasonable to expect students' pipelines to produce it as an alternative to the built-in module.
- Where exactly to set the three paralog-safety thresholds (score/coverage, score margin, minimum number of agreeing representatives) — a starting point from Phase 1's offline calibration, but worth revisiting once tested against real student assemblies.
- How aggressively to pre-cluster the marker-gene reference set before shipping it, trading index size/build time against how much representative diversity is kept for the multi-representative agreement check.
- What binary format the shipped seed index uses, and whether the build script's output should be versioned/checked into the repo directly or regenerated by a CI step, given it only needs rebuilding when the raw reference file or clustering/index parameters change.
- Which NCBI taxonomy dump version to resolve reference taxIDs against at build time, and how to handle any reference taxIDs that no longer resolve cleanly (merged/deleted/renamed nodes) if the dump postdates the reference file.
- What rank the marker-gene consensus lineage needs to resolve to before a bin is flagged as a chimerism candidate, versus treated as normal for a taxonomically well-resolved marker set — a starting point from Phase 1's calibration work, not a fixed rule.
- How many binning tools' worth of contig→bin tables the cross-tool reconciliation view needs to handle cleanly before the comparison UI (side-by-side views, colour coding) becomes unreadable, and what the fallback is beyond that (e.g. a summary matrix instead of full side-by-side).
- Whether the `beforeunload` guard should also cover filter/search state, or only actual contig reassignment, given the former is trivially reproducible and the latter isn't.
