#!/usr/bin/env node
'use strict';

// Step 1 of the offline marker-gene build pipeline (see
// docs/phase1-investigation.md §4). Pre-clusters reference-data/scg40_raw.fasta
// within each COG family at ~90% identity via a fast approximate method
// (minhash / k-mer-profile similarity — exact all-pairs alignment across
// ~1,700 sequences/family x 40 families is unnecessary for this cutoff).
//
// Input:  reference-data/scg40_raw.fasta
// Output: build/intermediate/scg40_clustered.fasta
//         build/intermediate/scg40_cluster_report.txt (representative
//         count per family, for eyeballing how much diversity survives)
//
// Not implemented yet — stub only, run order: 01 -> 02 -> 03.

throw new Error('build/01-cluster.js: not implemented yet');
