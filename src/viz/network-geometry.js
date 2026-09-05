(function () {
  'use strict';

// Pure layout for the cross-tool reconciliation bipartite network (MAG hubs
// vs. contig leaves) — no DOM here, matching reconciliation-network.js's
// split of geometry from rendering. Three interchangeable algorithms
// (layoutNetwork dispatches by name): a MAG-centric ring (deterministic,
// cheap, the default), a petal layout (leaves clustered into their primary
// hub's own arc around a shared inner ring of hubs, better once there are
// enough secondary hubs that the centric ring's crossings get busy), and a
// force-directed spring embedder (for genuinely dense/tangled networks, at
// the cost of an iterative computation and losing the ring's fixed-angle
// predictability).

/**
 * MAG-centric ring layout — the default arrangement. app.js's
 * buildMagNeighborhood always designates one hub as "central" (the MAG the
 * student selected, via opts.centralHubId); this layout puts that hub
 * literally at the canvas centre, every other ("secondary") hub it's
 * contending with on an outer ring, and contigs on a middle ring between
 * the two — clustered into a small angular sector near whichever secondary
 * hub they're contended with (same clustering idea as layoutPetal below),
 * so edges fan out locally instead of crossing the whole diagram.
 *
 * A leaf with only one hub (uncontended) doesn't belong on the "between"
 * ring at all, since there's nothing for it to sit between: one core to
 * the central hub sits close to the centre instead; one core to a
 * secondary hub (only possible when app.js's "show contigs for connected
 * MAGs" toggle pulls in that MAG's own uncontended contigs) sits in a
 * small arc just inside that hub's position on the outer ring.
 *
 * @param {string[]} hubIds
 * @param {{id:string, hubIds:string[]}[]} leaves - hubIds = distinct MAGs this leaf has an edge to
 * @param {{width:number, height:number, centralHubId?:string}} opts -
 *   centralHubId falls back to the first (sorted) hub id if omitted or not
 *   among hubIds, so the function stays well-defined without one
 * @returns {{hubPositions: Map<string,{x:number,y:number,angle:number}>,
 *            leafPositions: Map<string,{x:number,y:number,angle:number}>,
 *            orderedLeaves: object[]}}
 */
function layoutCentricRing(hubIds, leaves, opts) {
  const { width, height } = opts;
  const cx = width / 2, cy = height / 2;
  const outerR = Math.max(1, Math.min(width, height) / 2 - Math.max(24, Math.min(width, height) * 0.08));
  const midR = outerR * 0.55; // "between" ring for contended leaves
  const centreLeafR = outerR * 0.18; // near-centre ring for leaves core to the central hub
  const nearHubR = outerR * 0.85; // just inside the outer ring, for leaves core to a secondary hub

  const sortedHubIds = [...hubIds].sort();
  const centralHubId = hubIds.includes(opts.centralHubId) ? opts.centralHubId : sortedHubIds[0];
  const secondaryHubIds = sortedHubIds.filter((id) => id !== centralHubId);

  const hubPositions = new Map();
  if (centralHubId !== undefined) hubPositions.set(centralHubId, { x: cx, y: cy, angle: 0 });
  const hubAngle = new Map();
  secondaryHubIds.forEach((id, i) => {
    const angle = (i / Math.max(1, secondaryHubIds.length)) * 2 * Math.PI - Math.PI / 2;
    hubAngle.set(id, angle);
    hubPositions.set(id, { x: cx + outerR * Math.cos(angle), y: cy + outerR * Math.sin(angle), angle });
  });

  const contendedLeaves = [];
  const centreOnlyLeaves = [];
  const secondaryOnlyByHub = new Map();
  for (const leaf of leaves) {
    if (leaf.hubIds.length > 1) contendedLeaves.push(leaf);
    else if (leaf.hubIds[0] === centralHubId) centreOnlyLeaves.push(leaf);
    else {
      const hubId = leaf.hubIds[0];
      if (!secondaryOnlyByHub.has(hubId)) secondaryOnlyByHub.set(hubId, []);
      secondaryOnlyByHub.get(hubId).push(leaf);
    }
  }

  const leafPositions = new Map();
  const orderedLeaves = [];

  centreOnlyLeaves.sort((a, b) => a.id.localeCompare(b.id));
  centreOnlyLeaves.forEach((leaf, i) => {
    const angle = (i / Math.max(1, centreOnlyLeaves.length)) * 2 * Math.PI - Math.PI / 2;
    leafPositions.set(leaf.id, { x: cx + centreLeafR * Math.cos(angle), y: cy + centreLeafR * Math.sin(angle), angle });
    orderedLeaves.push(leaf);
  });

  const nearHubSpan = secondaryHubIds.length ? (2 * Math.PI / secondaryHubIds.length) * 0.5 : 2 * Math.PI;
  for (const [hubId, group] of secondaryOnlyByHub) {
    const center = hubAngle.get(hubId) ?? 0;
    group.sort((a, b) => a.id.localeCompare(b.id));
    group.forEach((leaf, i) => {
      const t = group.length === 1 ? 0.5 : i / (group.length - 1);
      const angle = center - nearHubSpan / 2 + t * nearHubSpan;
      leafPositions.set(leaf.id, { x: cx + nearHubR * Math.cos(angle), y: cy + nearHubR * Math.sin(angle), angle });
      orderedLeaves.push(leaf);
    });
  }

  // Contended leaves cluster near whichever secondary hub they contend
  // with — the smallest-id secondary among their hubs, for a stable,
  // deterministic grouping when a leaf touches more than two hubs.
  const byPrimarySecondary = new Map();
  for (const leaf of contendedLeaves) {
    const secondaries = leaf.hubIds.filter((id) => id !== centralHubId).sort();
    const primary = secondaries[0] ?? sortedHubIds[0];
    if (!byPrimarySecondary.has(primary)) byPrimarySecondary.set(primary, []);
    byPrimarySecondary.get(primary).push(leaf);
  }
  const midSpan = secondaryHubIds.length ? (2 * Math.PI / secondaryHubIds.length) * 0.8 : 2 * Math.PI;
  for (const [hubId, group] of byPrimarySecondary) {
    const center = hubAngle.get(hubId) ?? 0;
    group.sort((a, b) => a.id.localeCompare(b.id));
    group.forEach((leaf, i) => {
      const t = group.length === 1 ? 0.5 : i / (group.length - 1);
      const angle = center - midSpan / 2 + t * midSpan;
      leafPositions.set(leaf.id, { x: cx + midR * Math.cos(angle), y: cy + midR * Math.sin(angle), angle });
      orderedLeaves.push(leaf);
    });
  }

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

  const seed = layoutCentricRing(hubIds, leaves, opts);
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
  ring: layoutCentricRing,
  petal: layoutPetal,
  force: layoutForceDirected,
};

/** @param {'ring'|'petal'|'force'} algorithm */
function layoutNetwork(algorithm, hubIds, leaves, opts) {
  const fn = LAYOUT_ALGORITHMS[algorithm] || layoutCentricRing;
  return fn(hubIds, leaves, opts);
}

const exportsObj = { layoutCentricRing, layoutPetal, layoutForceDirected, layoutNetwork };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.networkGeometry = exportsObj;
}
})();
