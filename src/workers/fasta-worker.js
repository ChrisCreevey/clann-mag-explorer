'use strict';

// Runs the streaming FASTA parse off the main thread (Phase 2 performance
// follow-up — see docs/phase1-investigation.md "Phase 2 findings"). This
// doesn't reduce the total work, it just stops a multi-second-to-minutes
// parse from freezing the tab: the main thread posts a File, this worker
// streams it and posts per-contig records + periodic progress back, so the
// UI can render a live progress bar and stay responsive throughout.
//
// Marker-gene search itself runs in a separate pool of workers
// (src/workers/marker-search-worker.js, dispatched by src/app.js) so it can
// run in parallel across contigs — this worker keeps record.frames intact
// (rather than deleting it) so the main thread can hand it off.
//
// If any bin table was loaded, app.js hands over the union of every
// contigId those tables reference (`referencedContigIds`) — contigs the
// FASTA has but no bin table mentions are dropped right here, before
// translation/marker-search ever runs on them, rather than parsed in full
// and filtered out later. streamFasta's own summary always counts every
// record in the raw file (it increments before calling back), so this
// worker tracks its own kept-contig totals and reports those instead.
//
// importScripts (not <script> tags) since Workers don't share the page's
// DOM/script context; relative to this file's own location, not the page's.
importScripts('../model/dna-codes.js', '../model/translate.js', '../model/contig-stats.js', '../model/fasta-index.js');

const { streamFasta } = self.ClannMAG.fastaIndex;

const PROGRESS_EVERY = 200; // contigs between progress pings, not every single one

self.onmessage = async (e) => {
  const { file, referencedContigIds } = e.data;
  let count = 0;
  let keptContigCount = 0;
  let keptTotalLength = 0;
  try {
    const summary = await streamFasta(file, (record) => {
      if (referencedContigIds && !referencedContigIds.has(record.id)) return; // not in any loaded bin table — drop it
      count++;
      keptContigCount++;
      keptTotalLength += record.length;
      self.postMessage({ type: 'contig', record });
      if (count % PROGRESS_EVERY === 0) {
        self.postMessage({ type: 'progress', contigsSoFar: count });
      }
    });
    self.postMessage({
      type: 'done',
      summary: {
        ...summary,
        contigCount: keptContigCount,
        totalLength: keptTotalLength,
        excludedCount: summary.contigCount - keptContigCount,
      },
    });
  } catch (err) {
    self.postMessage({ type: 'error', message: err.message });
  }
};
