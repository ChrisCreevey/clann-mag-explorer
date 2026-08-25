'use strict';

// Standalone benchmark for the marker-gene search performance work (see
// docs/phase1-investigation.md and the maxSeedsPerSegment change in
// src/model/marker-genes.js / src/model/reduced-alphabet.js). Not part of
// the regular test/run.js suite — this measures wall-clock time, which
// makes it flaky as a pass/fail test but useful to re-run by hand whenever
// the search path changes.
//
// Usage:
//   node test/bench-marker-genes.js [path/to.fasta[.gz]] [--limit N] [--uncapped-sample N] [--workers N]
//
// Defaults to _realdata/contigs.binned.fasta.gz (a real ~45K-contig
// assembly already in this repo) if present, else examples/mini-assembly.fasta.
//
// Reports three numbers:
//   1. Single-threaded, capped (current default maxSeedsPerSegment=150) search time over the full set.
//   2. Single-threaded, uncapped (old behavior) search time over a sample, extrapolated to the full set
//      (running uncapped over the full set is what used to take ~20 minutes — too slow to just do here).
//   3. The same capped search run across a worker_threads pool, to measure the parallelism win from
//      splitting marker search into its own worker pool (src/app.js's MARKER_POOL_SIZE).

const fs = require('fs');
const path = require('path');
const zlib = require('zlib');
const { Worker, isMainThread, parentPort, workerData } = require('worker_threads');

const { translateSixFrames } = require('../src/model/translate');
const { loadMarkerGeneAssets, searchContigForMarkers, parseAssets } = require('../src/model/marker-genes');

const REPO_ROOT = path.join(__dirname, '..');

function bufToArrayBuffer(buf) {
  return buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength);
}

function loadAssetsFromDisk() {
  const dataDir = path.join(REPO_ROOT, 'data');
  const indexBuf = bufToArrayBuffer(fs.readFileSync(path.join(dataDir, 'scg40-index.bin')));
  const refSeqsBuf = bufToArrayBuffer(fs.readFileSync(path.join(dataDir, 'scg40-refseqs.bin')));
  const familyNames = JSON.parse(fs.readFileSync(path.join(dataDir, 'scg40-families.json'), 'utf8'));
  let thresholds;
  try {
    thresholds = JSON.parse(fs.readFileSync(path.join(dataDir, 'scg40-thresholds.json'), 'utf8'));
  } catch { /* optional */ }
  return parseAssets({ indexBuf, refSeqsBuf, familyNames, thresholds });
}

/** Minimal FASTA reader — good enough for a benchmark, not a replacement for src/model/fasta-index.js. */
function parseFasta(text) {
  const records = [];
  let id = null, chunks = [];
  const flush = () => { if (id) records.push({ id, seq: chunks.join('') }); };
  for (const line of text.split('\n')) {
    if (line.startsWith('>')) {
      flush();
      id = line.slice(1).trim().split(/\s+/)[0];
      chunks = [];
    } else if (line.length) {
      chunks.push(line.trim());
    }
  }
  flush();
  return records;
}

function loadFastaText(filePath) {
  const buf = fs.readFileSync(filePath);
  const raw = filePath.endsWith('.gz') ? zlib.gunzipSync(buf) : buf;
  return raw.toString('utf8');
}

function parseArgs(argv) {
  const args = { fastaPath: null, limit: Infinity, uncappedSample: 3000, workers: require('os').cpus().length };
  const rest = [];
  for (let i = 0; i < argv.length; i++) {
    if (argv[i] === '--limit') args.limit = Number(argv[++i]);
    else if (argv[i] === '--uncapped-sample') args.uncappedSample = Number(argv[++i]);
    else if (argv[i] === '--workers') args.workers = Number(argv[++i]);
    else if (argv[i] === '--skip-uncapped') args.skipUncapped = true;
    else if (argv[i] === '--skip-single-capped') args.skipSingleCapped = true;
    else rest.push(argv[i]);
  }
  args.fastaPath = rest[0] || null;
  return args;
}

function defaultFastaPath() {
  const realData = path.join(REPO_ROOT, '_realdata', 'contigs.binned.fasta.gz');
  if (fs.existsSync(realData)) return realData;
  return path.join(REPO_ROOT, 'examples', 'mini-assembly.fasta');
}

function fmt(ms) { return (ms / 1000).toFixed(1) + 's'; }

async function runMain() {
  const args = parseArgs(process.argv.slice(2));
  const fastaPath = args.fastaPath || defaultFastaPath();

  console.log(`Loading FASTA: ${fastaPath}`);
  const text = loadFastaText(fastaPath);
  let records = parseFasta(text);
  if (Number.isFinite(args.limit)) records = records.slice(0, args.limit);
  console.log(`${records.length.toLocaleString()} contigs loaded, ${args.workers} worker threads available`);

  console.log('Loading marker-gene assets from data/…');
  const assets = loadAssetsFromDisk();

  console.log('Translating six frames for every contig…');
  const t0 = performance.now();
  const framesByContig = records.map((r) => translateSixFrames(r.seq));
  console.log(`Translation: ${fmt(performance.now() - t0)}`);

  // 1. Single-threaded, capped (current default) — full set.
  if (!args.skipSingleCapped) {
    const t1 = performance.now();
    for (const frames of framesByContig) searchContigForMarkers(frames, assets);
    const elapsed = performance.now() - t1;
    console.log(`\n[1] Single-threaded, CAPPED (maxSeedsPerSegment=150), full set (${records.length.toLocaleString()} contigs): ${fmt(elapsed)}` +
      ` (${(elapsed / records.length).toFixed(2)} ms/contig)`);
  } else {
    console.log('\n[1] Skipped (--skip-single-capped)');
  }

  // 2. Single-threaded, uncapped (old behavior) — sample only, extrapolated.
  if (!args.skipUncapped) {
    const sampleSize = Math.min(args.uncappedSample, framesByContig.length);
    const sample = framesByContig.slice(0, sampleSize);
    const t2 = performance.now();
    for (const frames of sample) searchContigForMarkers(frames, assets, { maxSeedsPerSegment: Infinity });
    const elapsed = performance.now() - t2;
    const perContig = elapsed / sampleSize;
    console.log(`\n[2] Single-threaded, UNCAPPED (old behavior), sample of ${sampleSize.toLocaleString()} contigs: ${fmt(elapsed)}` +
      ` (${perContig.toFixed(2)} ms/contig)` +
      ` -> extrapolated full set (${records.length.toLocaleString()} contigs): ~${fmt(perContig * records.length)}`);
  } else {
    console.log('\n[2] Skipped (--skip-uncapped)');
  }

  // 3. Capped search across a worker_threads pool (mirrors src/app.js's marker-search-worker pool).
  {
    const poolSize = Math.max(1, args.workers);
    const t3 = performance.now();
    const elapsed = await runPooled(framesByContig, poolSize);
    console.log(`\n[3] Pooled across ${poolSize} worker threads, CAPPED, full set: ${fmt(elapsed)}` +
      ` (${(elapsed / records.length).toFixed(2)} ms/contig effective)`);
  }

  console.log('\nDone.');
}

function runPooled(framesByContig, poolSize) {
  return new Promise((resolve) => {
    const workers = Array.from({ length: poolSize }, () => new Worker(__filename, { workerData: { role: 'worker' } }));
    let nextIndex = 0, completed = 0;
    const t0 = performance.now();

    function dispatch(w) {
      if (nextIndex >= framesByContig.length) return;
      const idx = nextIndex++;
      w.postMessage({ frames: framesByContig[idx] });
    }

    workers.forEach((w) => {
      w.on('message', () => {
        completed++;
        if (completed === framesByContig.length) {
          const elapsed = performance.now() - t0;
          workers.forEach((ww) => ww.terminate());
          resolve(elapsed);
        } else {
          dispatch(w);
        }
      });
      dispatch(w);
    });
  });
}

if (isMainThread) {
  runMain().catch((err) => { console.error(err); process.exit(1); });
} else if (workerData && workerData.role === 'worker') {
  const assets = loadAssetsFromDisk();
  parentPort.on('message', ({ frames }) => {
    searchContigForMarkers(frames, assets);
    parentPort.postMessage({ done: true });
  });
}
