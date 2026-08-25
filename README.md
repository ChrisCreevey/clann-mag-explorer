# Clann MAG Explorer

## USE THE TOOL HERE: https://chriscreevey.github.io/clann-mag-explorer/

A browser-only tool for exploring and manually refining prokaryotic MAG (metagenome-assembled genome) binning
results from a single assembly. Fifth tool in the [Clann suite](https://chriscreevey.github.io/), alongside
[Clann Tree Viewer](https://chriscreevey.github.io/clann-tree-viewer/),
[Clann BLAST Explorer](https://chriscreevey.github.io/clann-blast-explorer/),
[Clann Pangenome Explorer](https://chriscreevey.github.io/clann-pangenome-explorer/), and
[Clann eDNA Explorer](https://chriscreevey.github.io/clann-edna-explorer/).

**Status: Phase 7 complete (streaming FASTA parsing + per-contig stats + marker-gene identification + bin
loading/summaries + cross-tool reconciliation + outlier/disagreement flagging + interactive reassignment), under
active development.** See [`clann-mag-explorer-brief.md`](clann-mag-explorer-brief.md) for the full design brief
and [`docs/phase1-investigation.md`](docs/phase1-investigation.md) for the current investigation/planning notes.

## Architecture

- Vanilla HTML/CSS/JS, no build step, no runtime dependencies, GPL-2.0.
- `index.html` + `styles/main.css` + `src/parsers/`, `src/model/`, `src/viz/` (plain `<script>` includes),
  matching the layout and design tokens of the other Clann tools.
- The streaming FASTA parse and per-contig stats run in a Web Worker (`src/workers/fasta-worker.js`) rather than
  the main thread, so loading a large assembly doesn't freeze the page. Every `src/parsers/`/`src/model/` module
  attaches itself to `self` (not `window`), so the same unmodified files load via `<script>` on the main thread
  and via `importScripts()` inside the worker.
- `build/` holds an **offline, Node-based** build pipeline (not part of the zero-dependency deployed site) that
  turns `reference-data/scg40_raw.fasta` into the static marker-gene search assets shipped in `data/`
  (`node build/01-cluster.js && node build/02-index.js`, ~25s total; `node build/03-taxonomy.js` separately,
  ~20s, downloads the current NCBI taxdump). `data/*.bin`/`data/*.json` are committed — GitHub Pages serves them
  directly, no CI build step. See `docs/phase1-investigation.md` §4-5 and "Phase 3/6 findings" for what each
  script does, the binary formats, and the performance tuning behind the current parameters.
- Marker-gene search (seed-and-extend against a Murphy10-reduced-alphabet index, BLOSUM62 extension, three
  paralog-safety checks) runs inside the same Worker as parsing, reusing the six-frame translation
  `contig-stats.js` already computes rather than recomputing it.
- Outlier/disagreement flagging (`src/model/outliers.js`, `src/model/marker-taxonomy.js`) combines
  composition/coverage centroid distance, marker-gene unique/redundant contribution, marker-gene taxonomic
  consistency (LCA of provenance taxIDs against `data/scg40-lineage.json`), Kraken2 per-contig disagreement, and
  Phase 5's cross-tool agreement into one ranked per-contig view.
- Interactive reassignment (`src/model/working-assignment.js`, `src/viz/scatter*.js`) is a live-editable bin
  assignment layered on top of whatever was loaded: an SVG scatter plot with rectangular drag-select (a
  simplified lasso — see `scatter-geometry.js`'s header for why), move/merge/new-bin actions, and bin summaries
  that recalculate immediately from the in-session edits. Arms a `beforeunload` guard on the first reassignment,
  not before.
- `test/` is a minimal zero-dependency harness (`node test/run.js`), matching the sibling tools.

## Repository layout

```
index.html                    Site shell
styles/main.css                Design tokens + layout (shared visual language with the other Clann tools)
src/parsers/                   Input-format parsers (contig->bin tables, coverage tables, Kraken2 .breport)
src/model/                     Core data model (FASTA index, contig stats, six-frame translation, taxonomy tree,
                                bin reconciliation, marker-gene search)
src/workers/                   Web Workers (streaming FASTA parse off the main thread)
src/viz/                       Visualisations (added as later phases land)
src/app.js                     App shell wiring
build/                         Offline build pipeline for the marker-gene reference assets (Node, not shipped)
data/                          Static assets shipped to the browser, built by build/ (index, ref seqs, lineage table)
reference-data/                Raw marker-gene reference FASTA (input to build/, not shipped to the browser)
docs/                          Investigation/planning docs
examples/                      Small example files for manual/browser testing
test/                          Zero-dependency test harness
```

## Running the tests

```
node test/run.js
```

## Licence

GPL-2.0, developed by [CreeveyLab](https://www.creeveylab.org/).
