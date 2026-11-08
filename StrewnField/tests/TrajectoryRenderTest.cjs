const fs = require('node:fs'),
  path = require('node:path'),
  assert = require('node:assert/strict'),
  crypto = require('node:crypto');
const native = require('@napi-rs/canvas'),
  T = require('../output/web/trajectory.js'),
  R = require('../output/web/report.js'),
  P = require('../output/web/print.js'),
  G = require('../output/web/graticule.js');
const html = fs.readFileSync(process.argv[2], 'utf8'),
  data = JSON.parse(html.match(/<script id="field-data" type="application\/json">([\s\S]*?)<\/script>/)[1]);
const run = data.runs[0],
  marker = `<script id="${run.flight.element}" type="application/json">`,
  start = html.indexOf(marker) + marker.length;
let flight = JSON.parse(html.slice(start, html.indexOf('</script>', start)));
if (flight.encoding === 'base64-deflate-json')
  flight = JSON.parse(require('node:zlib').inflateSync(Buffer.from(flight.payload, 'base64')));
const values = R.floats(
  require('node:zlib').inflateSync(Buffer.from(flight.samples, 'base64')).toString('base64'),
  true,
);
const model = T.prepare(flight, values, data.trace_columns.length, 120, 3, run.entry),
  stops = [
    [68, 1, 84],
    [59, 82, 139],
    [33, 145, 140],
    [94, 201, 98],
    [253, 231, 37],
  ];
const palette = (v) => {
  v = Math.max(0, Math.min(1, v)) * 4;
  const i = Math.min(3, Math.floor(v));
  return stops[i].map((c, k) => Math.round(c + (stops[i + 1][k] - c) * (v - i)));
};
const output = path.resolve(process.argv[3]);
fs.mkdirSync(output, { recursive: true });
(async () => {
  const tiles = [];
  for (const tile of data.tiles) tiles.push({ ...tile, image: await native.loadImage(tile.url) });
  const rows = model.tracks.flatMap((t) => t.rows.filter((r) => r[0] >= model.firstBreakup));
  const fallBounds = [
    [Math.min(...rows.map((r) => r[5])), Math.min(...rows.map((r) => r[6]))],
    [Math.max(...rows.map((r) => r[5])), Math.max(...rows.map((r) => r[6]))],
  ];
  for (const [name, width, height, bounds] of [
    ['entry', 1200, 760, model.bounds],
    ['fall', 1200, 760, fallBounds],
    ['mobile', 360, 600, fallBounds],
  ]) {
    const hashes = new Set();
    for (const [index, time] of [model.start, model.firstBreakup + 0.5, 60, model.end].entries()) {
      const c = native.createCanvas(width, height),
        ctx = c.getContext('2d'),
        p = P.projection(bounds, width, height, 45);
      ctx.fillStyle = '#47616a';
      ctx.fillRect(0, 0, width, height);
      for (const tile of tiles) {
        const n = 2 ** tile.z,
          [x, y] = p.world(tile.x / n, tile.y / n),
          [right, bottom] = p.world((tile.x + 1) / n, (tile.y + 1) / n);
        ctx.drawImage(tile.image, x, y, right - x, bottom - y);
      }
      G.draw(ctx, width, height, p.coordinate, p.point);
      const states = T.draw(
        ctx,
        model,
        time,
        (lat, lon) => {
          const [x, y] = p.point(lat, lon);
          return { x, y };
        },
        palette,
      );
      ctx.font = '13px Arial';
      ctx.fillStyle = '#fff';
      for (const m of model.markers) {
        const [x, y] = p.point(m.row[5], m.row[6]);
        if (x < 25 || x > width - 90 || y < 40 || y > height - 40) continue;
        ctx.beginPath();
        ctx.arc(x, y, m.kind === 'height' ? 3 : 6, 0, 2 * Math.PI);
        ctx.fill();
        if (m.kind === 'height') ctx.fillText(m.label, x + 8, y - 8);
      }
      const pixels = ctx.getImageData(0, 0, width, height).data;
      hashes.add(crypto.createHash('sha256').update(pixels).digest('hex'));
      assert.ok(states.length <= 120);
      fs.writeFileSync(path.join(output, `${name}-${index}.png`), c.toBuffer('image/png'));
    }
    assert.equal(hashes.size, 4, `${name}: all sampled animation frames must visibly change`);
  }
  for (const [name, width, height] of [
    ['profile-desktop', 780, 540],
    ['profile-mobile', 336, 400],
  ]) {
    const hashes = new Set();
    for (const [index, time] of [model.start, model.firstBreakup + 0.5, 60, model.end].entries()) {
      const c = native.createCanvas(width, height),
        ctx = c.getContext('2d'),
        states = T.frame(model, time);
      const profile = T.drawProfile(ctx, model, time, width, height, palette, states);
      assert.equal(profile.particles, states);
      for (const state of states) {
        const p = profile.point(state.row);
        assert.ok(
          p.x > 50 && p.x < width - 15 && p.y > 35 && p.y < height - 50,
          'Particle must remain framed',
        );
      }
      const pixels = ctx.getImageData(0, 0, width, height).data;
      let ink = 0;
      for (let i = 0; i < pixels.length; i += 4)
        if (Math.min(pixels[i], pixels[i + 1], pixels[i + 2]) < 200) ink++;
      assert.ok(ink > 1000, 'Profile axes, tracks and particles must not be blank');
      hashes.add(crypto.createHash('sha256').update(pixels).digest('hex'));
      fs.writeFileSync(path.join(output, `${name}-${index}.png`), c.toBuffer('image/png'));
    }
    assert.equal(hashes.size, 4, `${name}: all sampled frames must visibly change`);
  }
  console.log(
    'Native canvas frames verified at desktop/mobile sizes; this is not a browser layout test. ' + output,
  );
})().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
