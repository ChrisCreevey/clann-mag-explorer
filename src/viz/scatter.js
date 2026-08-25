(function () {
  'use strict';

// SVG scatter plot with rectangular drag-select (brief §Interactive
// reassignment — see scatter-geometry.js for why a rectangle, not a
// freeform lasso). Pure DOM/event wiring on top of scatter-geometry.js's
// scales and hit-testing; no rendering logic lives here that isn't about
// building/updating SVG elements or handling mouse events.

const SVG_NS = 'http://www.w3.org/2000/svg';

/**
 * @param {HTMLElement} container - emptied and filled with the plot
 * @param {{id:string, x:number, y:number, colorKey:string}[]} dataPoints
 * @param {{width?:number, height?:number, onSelectionChange:(ids:string[])=>void}} options
 * @returns {{setDataPoints:(pts:object[])=>void, getSelection:()=>string[], clearSelection:()=>void}}
 */
function createScatterPlot(container, dataPoints, options) {
  const { projectPoints, pointsInRect } = self.ClannMAG.scatterGeometry;
  const width = options.width || 600;
  const height = options.height || 400;
  const margin = 24;

  container.innerHTML = '';
  const svg = document.createElementNS(SVG_NS, 'svg');
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  svg.setAttribute('width', '100%');
  svg.style.display = 'block';
  svg.style.cursor = 'crosshair';
  container.appendChild(svg);

  const pointsGroup = document.createElementNS(SVG_NS, 'g');
  svg.appendChild(pointsGroup);
  const selectionRect = document.createElementNS(SVG_NS, 'rect');
  selectionRect.setAttribute('class', 'scatter-selection-rect');
  selectionRect.style.display = 'none';
  svg.appendChild(selectionRect);

  let currentDataPoints = dataPoints;
  let projected = [];
  let selectedIds = new Set();
  let dragStart = null;

  // Stable colour assignment per distinct colorKey (bin/MAG id) via CSS
  // custom properties on each point, cycling through a fixed hue set
  // rather than a data-driven scale — bin identity has no natural order
  // to map onto a gradient.
  const HUES = [200, 20, 150, 280, 50, 320, 100, 250, 0, 170];
  const hueByKey = new Map();
  function hueFor(key) {
    if (!hueByKey.has(key)) hueByKey.set(key, HUES[hueByKey.size % HUES.length]);
    return hueByKey.get(key);
  }

  function render() {
    const { points } = projectPoints(currentDataPoints, { width, height, margin });
    projected = points.map((p, i) => ({ ...p, colorKey: currentDataPoints[i].colorKey }));

    pointsGroup.innerHTML = '';
    for (const p of projected) {
      const circle = document.createElementNS(SVG_NS, 'circle');
      circle.setAttribute('cx', p.x);
      circle.setAttribute('cy', p.y);
      circle.setAttribute('r', selectedIds.has(p.id) ? 4.5 : 3);
      circle.setAttribute('fill', `hsl(${hueFor(p.colorKey)}, 55%, ${selectedIds.has(p.id) ? '45%' : '60%'})`);
      circle.setAttribute('stroke', selectedIds.has(p.id) ? 'var(--ink)' : 'none');
      circle.setAttribute('stroke-width', '1.5');
      circle.dataset.id = p.id;
      pointsGroup.appendChild(circle);
    }
  }

  function svgPoint(evt) {
    const rect = svg.getBoundingClientRect();
    return {
      x: ((evt.clientX - rect.left) / rect.width) * width,
      y: ((evt.clientY - rect.top) / rect.height) * height,
    };
  }

  svg.addEventListener('mousedown', (evt) => {
    dragStart = svgPoint(evt);
    selectionRect.style.display = 'block';
  });
  svg.addEventListener('mousemove', (evt) => {
    if (!dragStart) return;
    const cur = svgPoint(evt);
    const x = Math.min(dragStart.x, cur.x), y = Math.min(dragStart.y, cur.y);
    selectionRect.setAttribute('x', x);
    selectionRect.setAttribute('y', y);
    selectionRect.setAttribute('width', Math.abs(cur.x - dragStart.x));
    selectionRect.setAttribute('height', Math.abs(cur.y - dragStart.y));
  });
  function endDrag(evt) {
    if (!dragStart) return;
    const cur = svgPoint(evt);
    selectedIds = new Set(pointsInRect(projected, { x0: dragStart.x, y0: dragStart.y, x1: cur.x, y1: cur.y }));
    dragStart = null;
    selectionRect.style.display = 'none';
    render();
    options.onSelectionChange([...selectedIds]);
  }
  svg.addEventListener('mouseup', endDrag);
  svg.addEventListener('mouseleave', (evt) => { if (dragStart) endDrag(evt); });

  render();

  return {
    setDataPoints(pts) { currentDataPoints = pts; render(); },
    getSelection() { return [...selectedIds]; },
    clearSelection() { selectedIds = new Set(); render(); options.onSelectionChange([]); },
  };
}

const exportsObj = { createScatterPlot };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.scatter = exportsObj;
}
})();
