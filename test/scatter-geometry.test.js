const { test, report, assert } = require('./harness');
const { linearScale, domainOf, pointsInRect, projectPoints } = require('../src/viz/scatter-geometry');

test('domainOf finds min/max, [0,0] for empty input', () => {
  assert.deepStrictEqual(domainOf([3, 1, 4, 1, 5]), [1, 5]);
  assert.deepStrictEqual(domainOf([]), [0, 0]);
});

test('linearScale maps the domain endpoints into (not exactly onto) the padded range', () => {
  const scale = linearScale(0, 10, 0, 100, 0.1);
  assert.ok(scale(0) > 0, 'padding should keep the minimum off the very edge');
  assert.ok(scale(10) < 100, 'padding should keep the maximum off the very edge');
  assert.ok(Math.abs(scale(5) - 50) < 1e-9, 'midpoint should still map to the midpoint');
});

test('linearScale does not divide by zero when every value is identical', () => {
  const scale = linearScale(5, 5, 0, 100);
  assert.ok(Number.isFinite(scale(5)));
});

test('pointsInRect selects only points inside the rectangle, any corner order', () => {
  const points = [{ id: 'a', x: 5, y: 5 }, { id: 'b', x: 50, y: 50 }, { id: 'c', x: 95, y: 95 }];
  const selected = pointsInRect(points, { x0: 0, y0: 0, x1: 60, y1: 60 });
  assert.deepStrictEqual(selected.sort(), ['a', 'b']);
});

test('pointsInRect normalizes a rectangle dragged in any direction (e.g. bottom-right to top-left)', () => {
  const points = [{ id: 'a', x: 30, y: 30 }];
  const selected = pointsInRect(points, { x0: 60, y0: 60, x1: 0, y1: 0 });
  assert.deepStrictEqual(selected, ['a']);
});

test('pointsInRect includes points exactly on the boundary', () => {
  const points = [{ id: 'a', x: 10, y: 10 }];
  const selected = pointsInRect(points, { x0: 0, y0: 0, x1: 10, y1: 10 });
  assert.deepStrictEqual(selected, ['a']);
});

test('projectPoints keeps ids attached and produces coordinates within the plot area', () => {
  const dataPoints = [{ id: 'a', x: 0, y: 0 }, { id: 'b', x: 10, y: 100 }];
  const { points } = projectPoints(dataPoints, { width: 400, height: 300, margin: 20 });
  assert.strictEqual(points.length, 2);
  assert.deepStrictEqual(points.map((p) => p.id), ['a', 'b']);
  for (const p of points) {
    assert.ok(p.x >= 0 && p.x <= 400);
    assert.ok(p.y >= 0 && p.y <= 300);
  }
});

test('projectPoints flips the y axis so a larger data value plots higher (smaller pixel y)', () => {
  const dataPoints = [{ id: 'low', x: 0, y: 0 }, { id: 'high', x: 0, y: 100 }];
  const { points } = projectPoints(dataPoints, { width: 400, height: 300, margin: 20 });
  const low = points.find((p) => p.id === 'low');
  const high = points.find((p) => p.id === 'high');
  assert.ok(high.y < low.y, 'higher data value should have a smaller pixel y (closer to the top)');
});

report();
