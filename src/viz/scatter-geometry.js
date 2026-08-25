(function () {
  'use strict';

// Pure coordinate-mapping and selection-hit-testing for the interactive
// reassignment scatter plot (brief §Interactive reassignment). Kept
// separate from scatter.js's actual SVG/DOM construction so the geometry
// — the part with real logic worth getting right — can be unit tested
// without a browser DOM; scatter.js is a thin rendering/event-wiring
// layer on top of these functions.
//
// Scope decision: a rectangular drag-select, not a freeform lasso. A true
// lasso (arbitrary polygon, point-in-polygon hit testing) is a real
// feature, not just a naming difference, and a rectangle covers the same
// "select a visual cluster" need for a teaching tool's scatter plot at a
// fraction of the implementation/testing cost. Documented here rather
// than silently simplified.

/**
 * Linear scale from a data domain to a pixel range, with a fractional
 * padding fraction so points don't sit flush against the plot edge.
 * @returns {(value:number) => number}
 */
function linearScale(domainMin, domainMax, rangeMin, rangeMax, paddingFraction = 0.05) {
  const span = domainMax - domainMin;
  const pad = span === 0 ? 1 : span * paddingFraction; // guard: a single distinct value shouldn't divide by zero
  const paddedMin = domainMin - pad, paddedMax = domainMax + pad;
  const paddedSpan = paddedMax - paddedMin || 1;
  return (value) => rangeMin + ((value - paddedMin) / paddedSpan) * (rangeMax - rangeMin);
}

/**
 * @param {number[]} values
 * @returns {[number, number]} [min, max], both 0 if values is empty
 */
function domainOf(values) {
  if (values.length === 0) return [0, 0];
  let min = values[0], max = values[0];
  for (const v of values) { if (v < min) min = v; if (v > max) max = v; }
  return [min, max];
}

/**
 * @param {{id:string, x:number, y:number}[]} points - already in pixel space
 * @param {{x0:number, y0:number, x1:number, y1:number}} rect - drag rectangle, any corner order
 * @returns {string[]} ids of points inside the (normalized) rectangle
 */
function pointsInRect(points, rect) {
  const xMin = Math.min(rect.x0, rect.x1), xMax = Math.max(rect.x0, rect.x1);
  const yMin = Math.min(rect.y0, rect.y1), yMax = Math.max(rect.y0, rect.y1);
  return points
    .filter((p) => p.x >= xMin && p.x <= xMax && p.y >= yMin && p.y <= yMax)
    .map((p) => p.id);
}

/**
 * Maps data points to pixel-space points ready for rendering/selection,
 * given a plotting area and margins.
 * @param {{id:string, x:number, y:number}[]} dataPoints
 * @param {{width:number, height:number, margin:number}} plotArea
 * @returns {{points:{id:string,x:number,y:number}[], xScale:Function, yScale:Function}}
 */
function projectPoints(dataPoints, { width, height, margin }) {
  const [xMin, xMax] = domainOf(dataPoints.map((p) => p.x));
  const [yMin, yMax] = domainOf(dataPoints.map((p) => p.y));
  const xScale = linearScale(xMin, xMax, margin, width - margin);
  // SVG y grows downward; flip so larger data values plot higher on screen.
  const yScale = linearScale(yMin, yMax, height - margin, margin);
  const points = dataPoints.map((p) => ({ id: p.id, x: xScale(p.x), y: yScale(p.y) }));
  return { points, xScale, yScale };
}

const exportsObj = { linearScale, domainOf, pointsInRect, projectPoints };
if (typeof module !== 'undefined' && module.exports) module.exports = exportsObj;
if (typeof self !== 'undefined') {
  self.ClannMAG = self.ClannMAG || {};
  self.ClannMAG.scatterGeometry = exportsObj;
}
})();
