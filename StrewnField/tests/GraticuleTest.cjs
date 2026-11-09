const assert = require('node:assert/strict');
const G = require('../output/web/graticule.js');
const P = require('../output/web/print.js');
assert.equal(G.ticks(0.7, 0.9, 1000, 140).step, 0.05);
assert.equal(G.ticks(0.7, 0.8, 1000, 140).step, 0.02);
assert.equal(G.ticks(0.7, 0.9, 350, 140).step, 0.1);
for (const args of [
  [0, 0, 300, 120],
  [NaN, 1, 300, 120],
  [0, 1, 0, 120],
  [0, Infinity, 500, 100],
])
  assert.deepEqual(G.ticks(...args).values, []);
for (const [v, step, axis, want] of [
  [49.8, 0.02, 'latitude', '49.80\u00b0 N'],
  [-21.2, 0.05, 'latitude', '21.20\u00b0 S'],
  [0.8, 0.02, 'longitude', '0.80\u00b0 E'],
  [-66.1, 0.1, 'longitude', '66.1\u00b0 W'],
  [181, 1, 'longitude', '179\u00b0 W'],
  [-181, 1, 'longitude', '179\u00b0 E'],
  [180, 1, 'longitude', '180\u00b0'],
  [-1e-14, 0.02, 'longitude', '0.00\u00b0'],
  [49.800125, 0.000005, 'latitude', '49.800125\u00b0 N'],
])
  assert.equal(G.label(v, step, axis), want);
const bounds = [
  [
    [49.78, 0.72],
    [49.86, 0.9],
  ],
  [
    [-21.4, 22.5],
    [-21, 24],
  ],
  [
    [84, -20],
    [85, 20],
  ],
  [
    [0, 179],
    [1, 181],
  ],
  [
    [-80, -180],
    [80, 180],
  ],
  [
    [49.8001, 0.8001],
    [49.8002, 0.8002],
  ],
];
for (const b of bounds)
  for (const [w, h] of [
    [1600, 820],
    [650, 370],
    [320, 500],
  ]) {
    const p = P.projection(b, w, h),
      g = G.layout(w, h, p.coordinate, p.point);
    for (const [axis, ticks] of [
      ['longitude', g.longitudes],
      ['latitude', g.latitudes],
    ]) {
      assert.ok(ticks.length > 0 && ticks.length < 64);
      assert.equal(
        new Set(ticks.map((t) => t.label)).size,
        ticks.length,
        'Labels must not collapse to the same rounded coordinate',
      );
      for (const t of ticks) {
        assert.ok(t.position >= 0 && t.position <= (axis === 'longitude' ? w : h));
        const back =
          axis === 'longitude' ? p.coordinate(t.position, h / 2).lng : p.coordinate(w / 2, t.position).lat;
        assert.ok(Math.abs(back - t.value) < 1e-9, 'Grid lines must use the map projection');
      }
    }
    const labels = [],
      ctx = {
        save() {},
        restore() {},
        beginPath() {},
        rect() {},
        clip() {},
        moveTo() {},
        lineTo() {},
        setLineDash() {},
        stroke() {},
        fillRect() {},
        measureText: (s) => ({ width: s.length * 7 }),
        fillText: (text, x, y) => labels.push({ text, x, y, width: text.length * 7 }),
      };
    G.draw(ctx, w, h, p.coordinate, p.point);
    assert.ok(labels.length >= 2, 'Both coordinate axes must have readable labels');
    for (const t of labels) assert.ok(t.x > 0 && t.x + t.width < w && t.y > 0 && t.y < h);
    for (let i = 0; i < labels.length; i++)
      for (let j = i + 1; j < labels.length; j++) {
        const a = labels[i],
          c = labels[j];
        assert.ok(
          Math.abs(a.y - c.y) > 20 || a.x + a.width + 8 < c.x || c.x + c.width + 8 < a.x,
          'Coordinate labels must not overlap',
        );
      }
  }
console.log('Graticule spacing, hemispheres, dateline, projection and desktop/mobile label tests passed.');
