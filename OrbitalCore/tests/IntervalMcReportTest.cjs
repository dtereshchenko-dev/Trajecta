const { test } = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const R = require('../src/export/web/interval-mc.js');
test('print styles keep the heading legible and omit map action cells', () => {
  const web = path.join(__dirname, '../src/export/web');
  const shared = fs.readFileSync(path.join(web, 'map-print.css'), 'utf8');
  const mc = fs.readFileSync(path.join(web, 'interval-mc.css'), 'utf8');
  const common = fs.readFileSync(path.join(__dirname, '../../StrewnField/output/web/report-common.css'), 'utf8');
  assert.match(shared, /body\s*>\s*header\s*\{[^}]*background:\s*white\s*;/);
  assert.match(shared, /body\s+\.table-wrap\s*\{[^}]*max-height:\s*none\s*;/);
  assert.match(mc, /#entry-table th:last-child,\s*#entry-table td:last-child\s*\{\s*display:\s*none\s*;/);
  assert.match(common, /details:not\(\[open\]\)\s*\{\s*display:\s*none\s*;/);
  assert.match(common, /details\[open\]\s*>\s*summary\s*\{\s*display:\s*block\s*;/);
  assert.match(common, /details\[open\]\s+pre\s*\{[^}]*max-height:\s*none\s*;[^}]*overflow:\s*visible\s*;/);
});
test('UTC display preserves milliseconds and day rollover', () => {
  assert.equal(R.utc(54746 + (2 * 3600 + 45 * 60 + 29.23) / 86400), '2008-10-07 02:45:29.230');
  assert.equal(R.utc(54746 + 0.999999999), '2008-10-08 00:00:00.000');
  assert.equal(R.utc(null), 'UTC unavailable');
});
test('dateline corridor stays narrow without changing stored longitude', () => {
  const points = [{ longitude_deg: 179.9 }, { longitude_deg: -179.8 }];
  const centre = R.longitudeFrame(points),
    longitudes = points.map((p) => R.unwrap(p.longitude_deg, centre));
  assert.ok(Math.abs(longitudes[1] - longitudes[0]) < 0.31);
  assert.equal(points[1].longitude_deg, -179.8);
  assert.equal(R.longitudeFrame([]), 0);
});
test('event filter includes nominal separately and sorting does not mutate archive', () => {
  const points = [
    { clone: 7, event: 0, epoch: 3 },
    { clone: -1, event: 0, epoch: 2 },
    { clone: 2, event: 1, epoch: 1 },
  ];
  assert.deepEqual(
    R.ordered(points, '0', 'clone').map((p) => p.clone),
    [-1, 7],
  );
  assert.deepEqual(
    R.ordered(points, 'all', 'time').map((p) => p.clone),
    [2, -1, 7],
  );
  assert.equal(points[0].clone, 7);
  assert.deepEqual(R.ordered(points, '2', 'time'), []);
});
test('nominal markers are drawn last without changing the source order', () => {
  const points = [{ clone: -1, epoch: 1 }, { clone: 5, epoch: 2 }, { clone: -1, epoch: 3 }, { clone: 2, epoch: 4 }];
  const drawn = R.nominalLast(points);
  assert.deepEqual(drawn.map((p) => p.clone), [5, 2, -1, -1]);
  assert.deepEqual(points.map((p) => p.clone), [-1, 5, -1, 2]);
  assert.equal(drawn[2], points[0]);
  assert.deepEqual(R.nominalLast(points.filter((p) => p.clone >= 0)).map((p) => p.clone), [5, 2]);
});
test('clone and nominal visibility are independent and respect event selection', () => {
  const points = [
    { clone: -1, event: 0, epoch: 2 },
    { clone: 0, event: 0, epoch: 1 },
    { clone: 1, event: 1, epoch: 3 },
  ];
  assert.deepEqual(R.mapPoints(points, 'all', false, true).map((p) => p.clone), [-1]);
  assert.deepEqual(R.mapPoints(points, 'all', true, false).map((p) => p.clone), [0, 1]);
  assert.deepEqual(R.mapPoints(points, '0', true, true).map((p) => p.clone), [0, -1]);
  assert.deepEqual(R.mapPoints(points, 'all', false, false), []);
});
test('semilinear curves close only when all boundary crossings exist', () => {
  const p = (latitude_deg, longitude_deg) => ({ latitude_deg, longitude_deg });
  const closed = R.contourSegments([p(0, 179.9), p(1, -179.9), p(2, -179.8)], 180);
  assert.equal(closed.length, 1);
  assert.equal(closed[0].length, 4);
  assert.deepEqual(closed[0][0], closed[0][3]);
  assert.ok(Math.abs(closed[0][1][1] - closed[0][0][1]) < 1);
  const open = R.contourSegments([p(0, 0), p(1, 1), null, p(2, 2), p(3, 3)], 0);
  assert.equal(open.length, 1);
  assert.deepEqual(open[0], [[2, 2], [3, 3], [0, 0], [1, 1]]);
  assert.deepEqual(R.contourSegments([null, p(0, 0), null], 0), []);
});
test('contour bounds use available crossings from complete contours only', () => {
  const a = { latitude_deg: 1, longitude_deg: 2 };
  const b = { latitude_deg: 3, longitude_deg: 4 };
  const contours = [{ points: [a, null, b] }];
  assert.deepEqual(R.contourPoints({ status: 'complete', contours }), [a, b]);
  assert.deepEqual(R.contourPoints({ status: 'unavailable', contours }), []);
  assert.deepEqual(R.contourPoints(null), []);
});
test('semilinear status does not claim old or partial results are probabilities', () => {
  assert.match(R.contourNote(null), /not calculated/);
  assert.match(R.contourNote({ status: 'unavailable', reason: 'no covariance' }), /no covariance/);
  assert.match(R.contourNote({ status: 'complete', contours: [{ sigma: 5, missed: 2, failed: 1 }] }), /5σ: 2 missed, 1 failed/);
});
test('target selects a matching basemap including polar and empty planetary reports', () => {
  assert.equal(R.mapProfile('Earth', []).kind, 'xyz');
  assert.equal(R.mapProfile('Earth', [{ latitude_deg: 86 }]).kind, 'grid');
  assert.match(R.mapProfile('Moon', []).url, /Moon\/EQ\/LRO_WAC_Mosaic_Global/);
  assert.equal(R.mapProfile('Moon', [{ latitude_deg: 89.9 }]).kind, 'wmts');
  assert.match(R.mapProfile('Mars', []).url, /Mars\/EQ\/Mars_Viking_MDIM21/);
  assert.equal(R.mapProfile('Neptune', []).kind, 'grid');
  assert.equal(R.mapProfile('unknown', []).kind, 'grid');
  assert.equal(R.mapProfile('Moon', []).radius, 1737400);
  assert.equal(R.mapProfile('Mars', []).radius, 3396190);
});
test('planetary layers use the verified 2:1 WMTS pyramid with a retained overview', () => {
  const calls = [],
    L = {
      CRS: { EPSG4326: { code: 'EPSG:4326' } },
      tileLayer: (url, options) => {
        calls.push({ kind: 'xyz', url, options });
        return calls.at(-1);
      },
    };
  for (const target of ['Moon', 'Mars']) {
    const points = [{ latitude_deg: 45.123, longitude_deg: 179.8 }],
      copy = JSON.stringify(points);
    const profile = R.mapProfile(target, points),
      layer = R.imageryLayer(L, profile);
    assert.equal(layer.kind, 'xyz');
    assert.deepEqual(layer.options.bounds, [
      [-90, -180],
      [90, 180],
    ]);
    assert.match(layer.url, /^https:\/\/trek\.nasa\.gov\/tiles\//);
    assert.match(layer.url, /\/default028mm\/\{z\}\/\{y\}\/\{x\}\.jpg$/);
    assert.equal(layer.options.maxNativeZoom, target === 'Moon' ? 8 : 7);
    assert.equal(layer.options.updateWhenZooming, false);
    const overview = R.imageryLayer(L, profile, true);
    assert.equal(overview.url, layer.url);
    assert.equal(overview.options.minNativeZoom, 0);
    assert.equal(overview.options.maxNativeZoom, 0);
    assert.ok(overview.options.zIndex < layer.options.zIndex);
    assert.equal(JSON.stringify(points), copy, 'Report must preserve archived coordinates');
  }
  assert.equal(R.imageryLayer(L, R.mapProfile('Earth', [])).kind, 'xyz');
  assert.equal(R.imageryLayer(L, R.mapProfile('unknown', [])), null);
  assert.equal(calls.length, 5);
});
