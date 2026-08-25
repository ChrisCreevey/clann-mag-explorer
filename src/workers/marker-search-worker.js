'use strict';

// Runs marker-gene search (src/model/marker-genes.js) for one contig per
// message. Split out of fasta-worker.js so the main thread can run a pool
// of these (one per available core) and search contigs in parallel —
// fasta-worker.js stays single-threaded (streaming parse can't easily be
// sharded) but this, the actual bottleneck per contig, is embarrassingly
// parallel across contigs.
//
// importScripts relative to this file's own location, not the page's.
importScripts('../model/reduced-alphabet.js', '../model/blosum62.js', '../model/marker-genes.js');

const { loadMarkerGeneAssets, searchContigForMarkers } = self.ClannMAG.markerGenes;

const DATA_URL = new URL('../../data/', self.location).href;

const markerAssetsPromise = loadMarkerGeneAssets(DATA_URL).catch((err) => {
  self.postMessage({ type: 'marker-assets-error', message: err.message });
  return null;
});

self.onmessage = async (e) => {
  const { id, frames } = e.data;
  const assets = await markerAssetsPromise;
  const markerHits = (assets && frames) ? searchContigForMarkers(frames, assets) : [];
  self.postMessage({ id, markerHits });
};
