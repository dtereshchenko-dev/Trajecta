const { test } = require('node:test');
const assert = require('node:assert/strict');
const R = require('../src/export/web/interval-mc.js');
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
