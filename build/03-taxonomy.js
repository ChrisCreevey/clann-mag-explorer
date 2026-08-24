#!/usr/bin/env node
'use strict';

// Step 3 of the offline marker-gene build pipeline (see
// docs/phase1-investigation.md §4). Downloads the current NCBI taxdump
// (taxdump.tar.gz, build-time only, never checked in), extracts the
// distinct taxIDs present in reference-data/scg40_raw.fasta headers,
// resolves each against nodes.dmp/names.dmp, and handles merged/deleted
// taxIDs via merged.dmp/delnodes.dmp (map merged forward, warn-and-drop
// truly dead ones).
//
// Input:  reference-data/scg40_raw.fasta (for the taxID list)
//         NCBI taxdump, downloaded at build time
// Output: data/scg40-lineage.bin (or .json if small enough — TBD during
//         calibration; only covers the few thousand taxa actually
//         referenced, not the full taxdump)
//
// Not implemented yet — stub only, run order: 01 -> 02 -> 03.

throw new Error('build/03-taxonomy.js: not implemented yet');
