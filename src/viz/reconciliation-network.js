(function () {
  'use strict';

// Cross-tool reconciliation bipartite network: MAG hubs as central nodes,
// disputed contigs as leaves around them, one coloured edge per tool that
// voted a given contig into a given MAG (brief's core feature, "for each
// contig, how many of the loaded binning tools agree on where it
// belongs" — this makes that agreement structure visible directly, rather
// than read off an agreement-fraction table column). Pure DOM/SVG wiring
// on top of network-geometry.js's layouts, same split as scatter.js over
// scatter-geometry.js.
//
// Draggable, asymmetrically: dragging a MAG hub is "move this genome and
// everything currently blamed on it", so its connected leaves are dragged
// along with it as a rigid group — dragging a contig leaf is "just get
// this one contig out of the way", so it moves alone. That asymmetry is
// the point, not an inconsistency: a hub carries its leaves because the
// leaves' positions are being read *relative to* their hub, while a leaf
// has no dependents of its own to carry.

const SVG_NS = 'http://www.w3.org/2000/svg';
const HUES = [200, 20, 150, 280, 50, 320, 100, 250, 0, 170];

/**
 * @param {HTMLElement} container
 * @param {{hubs:{id:string,label:string}[], leaves:{id:string,hubIds:string[]}[],
 *           edges:{leafId:string,hubId:string,tool:string}[]}} data
 * @param {{width?:number, height?:number, algorithm?: 'ring'|'petal'|'force'}} options
 */
function createReconciliationNetwork(container, data, options = {}) {
  const { layoutNetwork } = self.ClannMAG.networkGeometry;
  const width = options.width || 640;
  const height = options.height || 640;
  const algorithm = options.algorithm || 'ring';

  container.innerHTML = '';

  const { hubs, leaves, edges } = data;
  if (leaves.length === 0) {
    container.innerHTML = '<div class="hint">No disputed contigs to show for the current filters.</div>';
    return;
  }

  const tools = [...new Set(edges.map((e) => e.tool))].sort();
  const hueByTool = new Map(tools.map((t, i) => [t, HUES[i % HUES.length]]));

  const wrap = document.createElement('div');
  wrap.className = 'network-wrap';
  const svg = document.createElementNS(SVG_NS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.setAttribute('width', '100%');
  svg.style.display = 'block';
  wrap.appendChild(svg);

  const legend = document.createElement('div');
  legend.className = 'network-legend';
  legend.innerHTML = `<span class="hint">Drag a MAG to move it with its contigs; drag a contig to move just that one.</span>` + tools
    .map((t) => `<span class="network-legend-item"><span class="network-swatch" style="background:hsl(${hueByTool.get(t)},60%,50%)"></span>${t}</span>`)
    .join('');
  wrap.appendChild(legend);
  container.appendChild(wrap);

  const { hubPositions, leafPositions } = layoutNetwork(algorithm, hubs.map((h) => h.id), leaves, { width, height });

  const edgesGroup = document.createElementNS(SVG_NS, 'g');
  const hubsGroup = document.createElementNS(SVG_NS, 'g');
  const leavesGroup = document.createElementNS(SVG_NS, 'g');
  svg.append(edgesGroup, hubsGroup, leavesGroup);

  const edgesByLeaf = new Map();
  const edgesByHub = new Map();
  const leafHubIds = new Map(leaves.map((l) => [l.id, l.hubIds]));

  for (const edge of edges) {
    const hp = hubPositions.get(edge.hubId);
    const lp = leafPositions.get(edge.leafId);
    if (!hp || !lp) continue;
    const line = document.createElementNS(SVG_NS, 'line');
    line.setAttribute('x1', hp.x); line.setAttribute('y1', hp.y);
    line.setAttribute('x2', lp.x); line.setAttribute('y2', lp.y);
    line.setAttribute('stroke', `hsl(${hueByTool.get(edge.tool)}, 60%, 50%)`);
    line.setAttribute('stroke-width', '1.4');
    line.setAttribute('opacity', '0.45');
    line.dataset.leaf = edge.leafId;
    line.dataset.hub = edge.hubId;
    edgesGroup.appendChild(line);
    if (!edgesByLeaf.has(edge.leafId)) edgesByLeaf.set(edge.leafId, []);
    edgesByLeaf.get(edge.leafId).push(line);
    if (!edgesByHub.has(edge.hubId)) edgesByHub.set(edge.hubId, []);
    edgesByHub.get(edge.hubId).push(line);
  }

  function setHighlight(activeLines) {
    const active = new Set(activeLines || []);
    edgesGroup.querySelectorAll('line').forEach((line) => {
      line.setAttribute('opacity', active.size === 0 ? '0.45' : active.has(line) ? '0.9' : '0.08');
    });
  }

  function setHubPosition(hubId, x, y) {
    hubPositions.set(hubId, { x, y });
    const el = hubEls.get(hubId);
    if (el) {
      el.circle.setAttribute('cx', x); el.circle.setAttribute('cy', y);
      el.label.setAttribute('x', x); el.label.setAttribute('y', y - 11);
    }
    for (const line of edgesByHub.get(hubId) || []) { line.setAttribute('x1', x); line.setAttribute('y1', y); }
  }

  function setLeafPosition(leafId, x, y) {
    leafPositions.set(leafId, { x, y });
    const el = leafEls.get(leafId);
    if (el) { el.circle.setAttribute('cx', x); el.circle.setAttribute('cy', y); }
    for (const line of edgesByLeaf.get(leafId) || []) { line.setAttribute('x2', x); line.setAttribute('y2', y); }
  }

  function svgPoint(evt) {
    const rect = svg.getBoundingClientRect();
    return {
      x: ((evt.clientX - rect.left) / rect.width) * width,
      y: ((evt.clientY - rect.top) / rect.height) * height,
    };
  }

  /**
   * Shared drag-start wiring: `getMovers()` is called fresh at drag-start
   * and returns the id->{x,y} snapshot of every node position this drag
   * should carry (just the dragged node itself for a leaf; the hub plus
   * every currently-connected leaf for a hub) — resolved at that moment,
   * not once up front, so it reflects any earlier drags. `apply(id,x,y)`
   * writes a moved position back (setHubPosition or setLeafPosition,
   * called per-mover with its own start position offset by the same drag
   * delta — the "rigid group move" for a hub).
   */
  function wireDrag(g, getMovers, apply) {
    g.style.cursor = 'grab';
    g.addEventListener('mousedown', (evt) => {
      evt.preventDefault();
      evt.stopPropagation();
      const start = svgPoint(evt);
      const movers = getMovers();
      g.style.cursor = 'grabbing';
      function onMove(moveEvt) {
        const cur = svgPoint(moveEvt);
        const dx = cur.x - start.x, dy = cur.y - start.y;
        for (const [id, p0] of movers) apply(id, p0.x + dx, p0.y + dy);
      }
      function onUp() {
        g.style.cursor = 'grab';
        window.removeEventListener('mousemove', onMove);
        window.removeEventListener('mouseup', onUp);
      }
      window.addEventListener('mousemove', onMove);
      window.addEventListener('mouseup', onUp);
    });
  }

  const hubEls = new Map();
  for (const hub of hubs) {
    const p = hubPositions.get(hub.id);
    if (!p) continue;
    const g = document.createElementNS(SVG_NS, 'g');
    const circle = document.createElementNS(SVG_NS, 'circle');
    circle.setAttribute('cx', p.x); circle.setAttribute('cy', p.y); circle.setAttribute('r', '7');
    circle.setAttribute('class', 'network-hub');
    const label = document.createElementNS(SVG_NS, 'text');
    label.setAttribute('x', p.x); label.setAttribute('y', p.y - 11);
    label.setAttribute('text-anchor', 'middle');
    label.setAttribute('class', 'network-hub-label');
    label.textContent = hub.label;
    g.append(circle, label);
    g.addEventListener('mouseenter', () => setHighlight(edgesByHub.get(hub.id)));
    g.addEventListener('mouseleave', () => setHighlight(null));
    const title = document.createElementNS(SVG_NS, 'title');
    title.textContent = `${hub.label} — ${(edgesByHub.get(hub.id) || []).length} vote(s) from disputed contigs shown. Drag to move it with its contigs.`;
    g.appendChild(title);
    hubsGroup.appendChild(g);
    hubEls.set(hub.id, { g, circle, label });
  }

  const leafEls = new Map();
  for (const leaf of leaves) {
    const p = leafPositions.get(leaf.id);
    if (!p) continue;
    const g = document.createElementNS(SVG_NS, 'g');
    const circle = document.createElementNS(SVG_NS, 'circle');
    circle.setAttribute('cx', p.x); circle.setAttribute('cy', p.y); circle.setAttribute('r', '3');
    circle.setAttribute('class', 'network-leaf');
    g.appendChild(circle);
    g.addEventListener('mouseenter', () => setHighlight(edgesByLeaf.get(leaf.id)));
    g.addEventListener('mouseleave', () => setHighlight(null));
    const title = document.createElementNS(SVG_NS, 'title');
    title.textContent = `${leaf.id} — voted into ${leaf.hubIds.length} different MAG(s) across tools. Drag to move just this contig.`;
    g.appendChild(title);
    leavesGroup.appendChild(g);
    leafEls.set(leaf.id, { g, circle });
  }

  // Wire dragging after all elements exist, since a hub's drag needs to
  // look up its connected leaves' *current* positions at drag-start time
  // (post any earlier drags), not the positions computed at initial layout.
  for (const hub of hubs) {
    const el = hubEls.get(hub.id);
    if (!el) continue;
    wireDrag(el.g, () => {
      const movers = new Map([[hub.id, { ...hubPositions.get(hub.id) }]]);
      for (const [leafId, hubIds] of leafHubIds) {
        if (hubIds.includes(hub.id)) movers.set(leafId, { ...leafPositions.get(leafId) });
      }
      return movers;
    }, (id, x, y) => (id === hub.id ? setHubPosition(id, x, y) : setLeafPosition(id, x, y)));
  }
  for (const leaf of leaves) {
    const el = leafEls.get(leaf.id);
    if (!el) continue;
    wireDrag(el.g, () => new Map([[leaf.id, { ...leafPositions.get(leaf.id) }]]), setLeafPosition);
  }
}

const exportsObj = { createReconciliationNetwork };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.reconciliationNetwork = exportsObj;
}
})();
