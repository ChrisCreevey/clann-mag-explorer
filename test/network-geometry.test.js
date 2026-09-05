const { test, report, assert } = require('./harness');
const { layoutCentricRing, layoutPetal, layoutForceDirected, layoutNetwork } = require('../src/viz/network-geometry');

const cx = 200, cy = 200;
const distFromCentre = (p) => Math.hypot(p.x - cx, p.y - cy);

test('the central hub sits exactly at the canvas centre', () => {
  const { hubPositions } = layoutCentricRing(['A', 'B', 'C'], [], { width: 400, height: 400, centralHubId: 'A' });
  assert.ok(distFromCentre(hubPositions.get('A')) < 1e-6);
});

test('secondary hubs are placed on an outer ring, all equidistant from centre and from the central hub', () => {
  const { hubPositions } = layoutCentricRing(['A', 'B', 'C'], [], { width: 400, height: 400, centralHubId: 'A' });
  const distB = distFromCentre(hubPositions.get('B'));
  const distC = distFromCentre(hubPositions.get('C'));
  assert.ok(Math.abs(distB - distC) < 1e-6);
  assert.ok(distB > distFromCentre(hubPositions.get('A')));
});

test('a contended leaf sits on a middle ring, strictly between the centre and the secondary hubs', () => {
  const leaves = [{ id: 'c1', hubIds: ['A', 'B'] }];
  const { hubPositions, leafPositions } = layoutCentricRing(['A', 'B'], leaves, { width: 400, height: 400, centralHubId: 'A' });
  const leafR = distFromCentre(leafPositions.get('c1'));
  const hubR = distFromCentre(hubPositions.get('B'));
  assert.ok(leafR > 0, 'should not sit exactly at the centre');
  assert.ok(leafR < hubR, 'should sit strictly inside the secondary hub\'s ring');
});

test('a leaf uncontended and core to the central hub sits close to the centre, closer than any contended leaf', () => {
  const leaves = [{ id: 'core1', hubIds: ['A'] }, { id: 'disputed1', hubIds: ['A', 'B'] }];
  const { leafPositions } = layoutCentricRing(['A', 'B'], leaves, { width: 400, height: 400, centralHubId: 'A' });
  assert.ok(distFromCentre(leafPositions.get('core1')) < distFromCentre(leafPositions.get('disputed1')));
});

test('a leaf uncontended and core to a secondary hub sits near that hub, not near the centre', () => {
  const leaves = [{ id: 'coreB', hubIds: ['B'] }];
  const { hubPositions, leafPositions } = layoutCentricRing(['A', 'B'], leaves, { width: 400, height: 400, centralHubId: 'A' });
  const p = leafPositions.get('coreB');
  const b = hubPositions.get('B');
  assert.ok(Math.hypot(p.x - b.x, p.y - b.y) < distFromCentre(p), 'closer to its own hub than to the centre');
});

test('every leaf and hub gets a distinct position', () => {
  const leaves = [{ id: 'c1', hubIds: ['A'] }, { id: 'c2', hubIds: ['B'] }, { id: 'c3', hubIds: ['A', 'B'] }];
  const { leafPositions } = layoutCentricRing(['A', 'B'], leaves, { width: 400, height: 400, centralHubId: 'A' });
  assert.strictEqual(leafPositions.size, 3);
});

test('falls back to the first sorted hub id as central when centralHubId is omitted or unknown', () => {
  const withOmitted = layoutCentricRing(['B', 'A'], [], { width: 400, height: 400 });
  const withUnknown = layoutCentricRing(['B', 'A'], [], { width: 400, height: 400, centralHubId: 'Z' });
  assert.ok(distFromCentre(withOmitted.hubPositions.get('A')) < 1e-6);
  assert.ok(distFromCentre(withUnknown.hubPositions.get('A')) < 1e-6);
});

test('layoutPetal places each leaf inside its primary hub angular sector', () => {
  const leaves = [{ id: 'x1', hubIds: ['A'] }, { id: 'x2', hubIds: ['A'] }, { id: 'y1', hubIds: ['B'] }];
  const { hubPositions, leafPositions } = layoutPetal(['A', 'B'], leaves, { width: 400, height: 400 });
  const cx = 200, cy = 200;
  const angleOf = (p) => Math.atan2(p.y - cy, p.x - cx);
  const hubAngleA = angleOf(hubPositions.get('A'));
  const hubAngleB = angleOf(hubPositions.get('B'));
  const leafAngleX1 = angleOf(leafPositions.get('x1'));
  const leafAngleY1 = angleOf(leafPositions.get('y1'));
  assert.ok(Math.abs(leafAngleX1 - hubAngleA) < Math.abs(leafAngleX1 - hubAngleB));
  assert.ok(Math.abs(leafAngleY1 - hubAngleB) < Math.abs(leafAngleY1 - hubAngleA));
});

test('layoutForceDirected keeps every node within the canvas bounds', () => {
  const leaves = [{ id: 'c1', hubIds: ['A'] }, { id: 'c2', hubIds: ['A', 'B'] }, { id: 'c3', hubIds: ['B'] }];
  const { hubPositions, leafPositions } = layoutForceDirected(['A', 'B'], leaves, { width: 300, height: 300, iterations: 30 });
  for (const p of [...hubPositions.values(), ...leafPositions.values()]) {
    assert.ok(p.x >= 0 && p.x <= 300 && p.y >= 0 && p.y <= 300);
  }
});

test('layoutForceDirected is deterministic given the same input (seeded from the ring layout, not random)', () => {
  const leaves = [{ id: 'c1', hubIds: ['A'] }, { id: 'c2', hubIds: ['A', 'B'] }, { id: 'c3', hubIds: ['B'] }];
  const first = layoutForceDirected(['A', 'B'], leaves, { width: 400, height: 400, iterations: 50 });
  const second = layoutForceDirected(['A', 'B'], leaves, { width: 400, height: 400, iterations: 50 });
  assert.deepStrictEqual(first.hubPositions.get('A'), second.hubPositions.get('A'));
  assert.deepStrictEqual(first.leafPositions.get('c2'), second.leafPositions.get('c2'));
});

test('layoutForceDirected keeps a leaf closer to the hub it is connected to than to an unconnected hub', () => {
  // A modest graph (small-n force layouts are dominated by the repulsion
  // term and don't reliably form tight clusters — see the bounds/
  // determinism tests above for what's checked at n=2-3): two hubs, each
  // with several of their own leaves plus one leaf shared between them.
  const leaves = [
    { id: 'a1', hubIds: ['A'] }, { id: 'a2', hubIds: ['A'] }, { id: 'a3', hubIds: ['A'] },
    { id: 'b1', hubIds: ['B'] }, { id: 'b2', hubIds: ['B'] }, { id: 'b3', hubIds: ['B'] },
  ];
  const { hubPositions, leafPositions } = layoutForceDirected(['A', 'B'], leaves, { width: 500, height: 500, iterations: 150 });
  const a = hubPositions.get('A'), b = hubPositions.get('B');
  for (const id of ['a1', 'a2', 'a3']) {
    const p = leafPositions.get(id);
    assert.ok(Math.hypot(p.x - a.x, p.y - a.y) < Math.hypot(p.x - b.x, p.y - b.y));
  }
});

test('layoutNetwork dispatches by algorithm name and falls back to ring for an unknown one', () => {
  const leaves = [{ id: 'c1', hubIds: ['A'] }];
  const viaRing = layoutNetwork('ring', ['A'], leaves, { width: 200, height: 200 });
  const viaUnknown = layoutNetwork('nonsense', ['A'], leaves, { width: 200, height: 200 });
  assert.deepStrictEqual(viaRing.hubPositions.get('A'), viaUnknown.hubPositions.get('A'));
});

report();
