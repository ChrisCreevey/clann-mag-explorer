#!/usr/bin/env node
'use strict';

// Step 3 of the offline marker-gene build pipeline (see
// docs/phase1-investigation.md §4, and "Phase 3 findings" for why this
// was deferred until Phase 6 actually needed it — completeness/redundancy
// and per-contig family tags don't need lineage, but Phase 6's marker-gene
// taxonomic-consistency check does). Downloads the current NCBI taxdump
// (build-time only, never checked in), resolves every distinct taxID
// referenced in reference-data/scg40_raw.fasta headers to its full
// ancestor chain, and ships a compact taxid->{parent,rank,name} table —
// just enough for src/model/taxonomy-tree.js's existing generic LCA logic
// to run at runtime against real NCBI lineage instead of a Kraken2
// breport's indentation-derived tree (see docs/phase1-investigation.md §2
// for why that tree needed a different ingestion path).
//
// Input:  reference-data/scg40_raw.fasta (for the taxID list)
//         NCBI taxdump, downloaded fresh at build time
// Output: data/scg40-lineage.json — {taxid, parentTaxid, rank, name} for
//         every referenced taxID and all of its ancestors up to root.
//         JSON, not binary: a few thousand nodes is small enough that the
//         binary-format tradeoffs flagged in docs/phase1-investigation.md
//         §5 don't apply here the way they do to the seed index.

const fs = require('fs');
const path = require('path');
const { execSync } = require('child_process');
const { parseFasta, parseHeader } = require('./02-index');

const TAXDUMP_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz';
const REFERENCE_FASTA = path.join(__dirname, '..', 'reference-data', 'scg40_raw.fasta');
const WORK_DIR = path.join(__dirname, 'intermediate', 'taxdump');
const OUT_PATH = path.join(__dirname, '..', 'data', 'scg40-lineage.json');

function downloadAndExtractTaxdump() {
  fs.mkdirSync(WORK_DIR, { recursive: true });
  const archivePath = path.join(WORK_DIR, 'taxdump.tar.gz');
  const nodesPath = path.join(WORK_DIR, 'nodes.dmp');
  if (fs.existsSync(nodesPath)) {
    console.log(`Reusing already-extracted taxdump at ${WORK_DIR}`);
    return;
  }
  console.log(`Downloading ${TAXDUMP_URL} ...`);
  execSync(`curl -sL --fail -o "${archivePath}" "${TAXDUMP_URL}"`, { stdio: 'inherit' });
  console.log('Extracting nodes.dmp, names.dmp, merged.dmp, delnodes.dmp ...');
  execSync(`tar -xzf "${archivePath}" -C "${WORK_DIR}" nodes.dmp names.dmp merged.dmp delnodes.dmp`, { stdio: 'inherit' });
}

/** taxdump .dmp files: rows separated by "\t|\n", fields by "\t|\t". */
function parseDmpLines(text) {
  return text.split('\t|\n').filter((line) => line.length > 0).map((line) => line.split('\t|\t'));
}

function loadNodes() {
  const text = fs.readFileSync(path.join(WORK_DIR, 'nodes.dmp'), 'utf8');
  const nodes = new Map(); // taxid -> {parentTaxid, rank}
  for (const fields of parseDmpLines(text)) {
    const taxid = Number(fields[0]);
    const parentTaxid = Number(fields[1]);
    const rank = fields[2];
    nodes.set(taxid, { parentTaxid, rank });
  }
  return nodes;
}

function loadScientificNames() {
  const text = fs.readFileSync(path.join(WORK_DIR, 'names.dmp'), 'utf8');
  const names = new Map(); // taxid -> scientific name
  for (const fields of parseDmpLines(text)) {
    const taxid = Number(fields[0]);
    const name = fields[1];
    const nameClass = fields[3];
    if (nameClass === 'scientific name') names.set(taxid, name);
  }
  return names;
}

function loadMerged() {
  const text = fs.readFileSync(path.join(WORK_DIR, 'merged.dmp'), 'utf8');
  const merged = new Map(); // old taxid -> new taxid
  for (const fields of parseDmpLines(text)) merged.set(Number(fields[0]), Number(fields[1]));
  return merged;
}

function loadDeletedSet() {
  const text = fs.readFileSync(path.join(WORK_DIR, 'delnodes.dmp'), 'utf8');
  return new Set(parseDmpLines(text).map((fields) => Number(fields[0])));
}

/**
 * NCBI's rank strings ("superkingdom", "phylum", ...) don't match the
 * single-letter rank codes taxonomy-tree.js's canonicalRank() expects
 * (Kraken2 .breport convention, brief's other input source for this same
 * tree). Mapped here, at build time, so both ingestion paths produce
 * rank codes in the same shape the tree/LCA logic already understands.
 * Ranks with no standard Kraken-style single-letter code fall back to a
 * literal-string "rank code" (canonicalRank treats any unrecognised
 * string as its own rankLetter, sub 0) — imprecise for depth-comparison
 * purposes but harmless for LCA, which only walks parent pointers.
 */
const RANK_CODE = {
  'superkingdom': 'D', 'kingdom': 'K', 'phylum': 'P', 'class': 'C',
  'order': 'O', 'family': 'F', 'genus': 'G', 'species': 'S',
  'no rank': 'R', 'root': 'R',
};
function rankCode(ncbiRank) { return RANK_CODE[ncbiRank] || ncbiRank; }

function main() {
  downloadAndExtractTaxdump();

  console.log(`Reading ${REFERENCE_FASTA} for referenced taxIDs...`);
  const records = parseFasta(fs.readFileSync(REFERENCE_FASTA, 'utf8')).map((r) => ({ ...r, ...parseHeader(r.header) }));
  const referencedTaxIds = new Set(records.map((r) => r.taxId));
  console.log(`${referencedTaxIds.size} distinct taxIDs referenced across ${records.length} sequences.`);

  const nodes = loadNodes();
  const names = loadScientificNames();
  const merged = loadMerged();
  const deleted = loadDeletedSet();

  const resolved = new Map(); // taxid -> {taxid, parentTaxid, rank, name}
  let mergedCount = 0, deadCount = 0;

  function resolveAncestorChain(startTaxid) {
    let taxid = startTaxid;
    if (merged.has(taxid)) { taxid = merged.get(taxid); mergedCount++; }
    if (!nodes.has(taxid)) {
      if (deleted.has(taxid)) deadCount++;
      return; // warn-and-drop, per docs/phase1-investigation.md's planned handling
    }
    while (taxid !== undefined && !resolved.has(taxid)) {
      const node = nodes.get(taxid);
      if (!node) break;
      resolved.set(taxid, { taxid, parentTaxid: node.parentTaxid, rank: rankCode(node.rank), name: names.get(taxid) || `taxid:${taxid}` });
      if (taxid === node.parentTaxid) break; // root (taxid 1) is its own parent
      taxid = node.parentTaxid;
    }
  }

  for (const taxid of referencedTaxIds) resolveAncestorChain(taxid);

  console.log(`Resolved ${resolved.size} nodes total (referenced taxIDs + all ancestors up to root).`);
  console.log(`${mergedCount} referenced taxIDs were merged forward; ${deadCount} referenced taxIDs no longer resolve (deleted) and were dropped.`);

  fs.mkdirSync(path.dirname(OUT_PATH), { recursive: true });
  fs.writeFileSync(OUT_PATH, JSON.stringify([...resolved.values()]));
  console.log(`Wrote ${OUT_PATH}`);
}

if (require.main === module) main();

module.exports = { parseDmpLines, rankCode };
