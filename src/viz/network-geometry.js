(function () {
  'use strict';

// Pure layout for the cross-tool reconciliation bipartite network (MAG hubs
// vs. contig leaves) — no DOM here, matching reconciliation-network.js's
// split of geometry from rendering. Three interchangeable algorithms
// (layoutNetwork dispatches by name): a plain ring (deterministic, cheap,
// good default for smaller networks), a petal layout (leaves clustered
// into their primary hub's own arc, better once there are enough hubs
// that a shared ring's crossings get busy), and a force-directed spring
// embedder (for genuinely dense/tangled networks, at the cost of an
// iterative computation and losing the ring's fixed-angle predictability).

/**
 * @param {string[]} hubIds
 * @param {{id:string, hubIds:string[]}[]} leaves - hubIds = distinct MAGs this leaf has an edge to
 * @param {{width:number, height:number, hubRadiusFraction?:number}} opts
 * @returns {{hubPositions: Map<string,{x:number,y:number,angle:number}>,
 *            leafPositions: Map<string,{x:number,y:number,angle:number}>,
 *            orderedLeaves: object[]}}
 */
function layoutBipartiteNetwork(hubIds, leaves, opts) {
  const { width, height, hubRadiusFraction = 0.32 } = opts;
  const cx = width / 2, cy = height / 2;
  const outerR = Math.max(1, Math.min(width, height) / 2 - Math.max(24, Math.min(width, height) * 0.08));
  const innerR = outerR * hubRadiusFraction;

  // Cluster leaves that share a primary (first, lexicographically smallest)
  // hub next to each other around the ring, so edges from one hub's
  // disputed contigs fan out together rather than crossing the whole
  // circle at random.
  const orderedLeaves = [...leaves].sort((a, b) => {
    const ah = [...a.hubIds].sort()[0] || '';
    const bh = [...b.hubIds].sort()[0] || '';
    return ah.localeCompare(bh) || a.id.localeCompare(b.id);
  });

  const hubPositions = new Map();
  const sortedHubIds = [...hubIds].sort();
  sortedHubIds.forEach((id, i) => {
    const angle = (i / Math.max(1, sortedHubIds.length)) * 2 * Math.PI - Math.PI / 2;
    hubPositions.set(id, { x: cx + innerR * Math.cos(angle), y: cy + innerR * Math.sin(angle), angle });
  });

  const leafPositions = new Map();
  orderedLeaves.forEach((leaf, i) => {
    const angle = (i / Math.max(1, orderedLeaves.length)) * 2 * Math.PI - Math.PI / 2;
    leafPositions.set(leaf.id, { x: cx + outerR * Math.cos(angle), y: cy + outerR * Math.sin(angle), angle });
  });

  return { hubPositions, leafPositions, orderedLeaves };
}

/**
 * Petal layout: same inner ring of hubs, but each leaf sits in its primary
 * hub's own angular sector of the outer ring rather than one shared
 * continuous ring — so a hub's disputed contigs fan out in a compact arc
 * right in front of it ("petal") instead of occupying a stretch of a
 * global ring that a leaf's *secondary* hub's edge still has to cross to
 * reach. Better than the plain ring once there are enough hubs that a
 * leaf's petal and its secondary-hub edge would otherwise cross several
 * unrelated petals.
 */
function layoutPetal(hubIds, leaves, opts) {
  const { width, height, hubRadiusFraction = 0.28 } = opts;
  const cx = width / 2, cy = height / 2;
  const outerR = Math.max(1, Math.min(width, height) / 2 - Math.max(24, Math.min(width, height) * 0.08));
  const innerR = outerR * hubRadiusFraction;

  const sortedHubIds = [...hubIds].sort();
  const hubAngle = new Map();
  sortedHubIds.forEach((id, i) => hubAngle.set(id, (i / Math.max(1, sortedHubIds.length)) * 2 * Math.PI - Math.PI / 2));

  const hubPositions = new Map();
  for (const id of sortedHubIds) {
    const angle = hubAngle.get(id);
    hubPositions.set(id, { x: cx + innerR * Math.cos(angle), y: cy + innerR * Math.sin(angle), angle });
  }

  const byPrimaryHub = new Map();
  for (const leaf of leaves) {
    const primary = [...leaf.hubIds].sort()[0];
    if (!byPrimaryHub.has(primary)) byPrimaryHub.set(primary, []);
    byPrimaryHub.get(primary).push(leaf);
  }

  const sectorSpan = (2 * Math.PI / Math.max(1, sortedHubIds.length)) * 0.8; // gap between petals so they read as distinct
  const leafPositions = new Map();
  const orderedLeaves = [];
  for (const [hubId, group] of byPrimaryHub) {
    const center = hubAngle.get(hubId) ?? 0;
    group.sort((a, b) => a.id.localeCompare(b.id));
    group.forEach((leaf, i) => {
      const t = group.length === 1 ? 0.5 : i / (group.length - 1);
      const angle = center - sectorSpan / 2 + t * sectorSpan;
      leafPositions.set(leaf.id, { x: cx + outerR * Math.cos(angle), y: cy + outerR * Math.sin(angle), angle });
      orderedLeaves.push(leaf);
    });
  }

  return { hubPositions, leafPositions, orderedLeaves };
}

/**
 * Force-directed (Fruchterman-Reingold-style spring embedder) layout, for
 * dense/messy networks where the fixed ring/petal layouts still cross a
 * lot of edges. Nodes repel each other, edges act like springs pulling
 * connected hub/leaf pairs together, and a mild pull-to-centre keeps the
 * whole graph from drifting off-canvas. Seeded from the ring layout
 * (rather than random positions) so it converges faster and produces a
 * broadly similar, only-locally-rearranged result on repeat runs, rather
 * than a different layout every time. Runs its iterations synchronously
 * and returns the final positions — not animated, since a one-off
 * "arrange" action doesn't need to be watched settling.
 */
function layoutForceDirected(hubIds, leaves, opts) {
  const { width, height, iterations = 150 } = opts;
  const cx = width / 2, cy = height / 2;

  const seed = layoutBipartiteNetwork(hubIds, leaves, opts);
  const nodes = new Map();
  for (const [id, p] of seed.hubPositions) nodes.set(id, { x: p.x, y: p.y });
  for (const [id, p] of seed.leafPositions) nodes.set(id, { x: p.x, y: p.y });

  const edges = [];
  for (const leaf of leaves) for (const hubId of leaf.hubIds) edges.push([leaf.id, hubId]);

  const ids = [...nodes.keys()];
  const k = Math.sqrt((width * height) / Math.max(1, ids.length)) * 0.9;
  let temperature = Math.min(width, height) / 10;

  for (let iter = 0; iter < iterations; iter++) {
    const disp = new Map(ids.map((id) => [id, { x: 0, y: 0 }]));

    for (let i = 0; i < ids.length; i++) {
      for (let j = i + 1; j < ids.length; j++) {
        const a = nodes.get(ids[i]), b = nodes.get(ids[j]);
        let dx = a.x - b.x, dy = a.y - b.y;
        const dist = Math.hypot(dx, dy) || 0.01;
        const force = (k * k) / dist;
        dx = (dx / dist) * force; dy = (dy / dist) * force;
        disp.get(ids[i]).x += dx; disp.get(ids[i]).y += dy;
        disp.get(ids[j]).x -= dx; disp.get(ids[j]).y -= dy;
      }
    }
    for (const [leafId, hubId] of edges) {
      const a = nodes.get(leafId), b = nodes.get(hubId);
      if (!a || !b) continue;
      let dx = a.x - b.x, dy = a.y - b.y;
      const dist = Math.hypot(dx, dy) || 0.01;
      const force = (dist * dist) / k;
      dx = (dx / dist) * force; dy = (dy / dist) * force;
      disp.get(leafId).x -= dx; disp.get(leafId).y -= dy;
      disp.get(hubId).x += dx; disp.get(hubId).y += dy;
    }

    for (const id of ids) {
      const n = nodes.get(id);
      const d = disp.get(id);
      const dlen = Math.hypot(d.x, d.y) || 0.01;
      n.x += (d.x / dlen) * Math.min(dlen, temperature);
      n.y += (d.y / dlen) * Math.min(dlen, temperature);
      n.x += (cx - n.x) * 0.01;
      n.y += (cy - n.y) * 0.01;
      n.x = Math.min(width - 10, Math.max(10, n.x));
      n.y = Math.min(height - 10, Math.max(10, n.y));
    }
    temperature *= 0.97;
  }

  const hubPositions = new Map(hubIds.map((id) => [id, { x: nodes.get(id).x, y: nodes.get(id).y }]));
  const leafPositions = new Map(leaves.map((l) => [l.id, { x: nodes.get(l.id).x, y: nodes.get(l.id).y }]));
  return { hubPositions, leafPositions, orderedLeaves: leaves };
}

const LAYOUT_ALGORITHMS = {
  ring: layoutBipartiteNetwork,
  petal: layoutPetal,
  force: layoutForceDirected,
};

/** @param {'ring'|'petal'|'force'} algorithm */
function layoutNetwork(algorithm, hubIds, leaves, opts) {
  const fn = LAYOUT_ALGORITHMS[algorithm] || layoutBipartiteNetwork;
  return fn(hubIds, leaves, opts);
}

const exportsObj = { layoutBipartiteNetwork, layoutPetal, layoutForceDirected, layoutNetwork };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.networkGeometry = exportsObj;
}
})();
