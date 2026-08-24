#!/usr/bin/env node
'use strict';

// Step 2 of the offline marker-gene build pipeline (see
// docs/phase1-investigation.md §4-5). Builds the seed index from the
// clustered reference set: reduced-alphabet translation (Murphy-10 by
// default), k=5 windowing (k=6 fallback if calibration shows too many
// collisions), hashed to integer keys, CSR-style typed arrays.
//
// Input:  build/intermediate/scg40_clustered.fasta
// Output: data/scg40-index.bin
//         data/scg40-refseqs.bin (real residues, for extension scoring)
//
// Not implemented yet — stub only, run order: 01 -> 02 -> 03.

throw new Error('build/02-index.js: not implemented yet');
