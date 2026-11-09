const assert = require('node:assert/strict');
const T = require('../output/web/trajectory.js');
const rows = (time, height, mass, lon, speed = 10) => [time, height, mass, speed, 0, 49, lon, speed];
const flight = { verified: true, traces: [], breakups: [] },
  values = [];
function add(id, parent, generation, end, a, b) {
  const t = {
    id,
    parent_id: parent,
    generation,
    end_event: end,
    multiplicity: 1,
    count: 2,
    offset: values.length,
    first: a,
    last: b,
  };
  flight.traces.push(t);
  values.push(...a, ...b);
  return t;
}
add(0, 0, 0, 'breakup', rows(0, 100, 100, 179), rows(5, 30, 90, 179.5));
flight.breakups.push({ fragment_id: 0, state: rows(5, 30, 90, 179.5), pressure_pa: 4e6 });
add(1, 0, 1, 'ground', rows(5, 30, 30, 179.5), rows(15, 0, 5, -179.5, 0.1));
add(2, 0, 1, 'mass_cutoff', rows(5, 30, 1, 179.5), rows(7, 20, 0.001, 179.8, 1));
const prepare = (limit = 120) => T.prepare(flight, Float32Array.from(values), 8, limit);
const model = prepare();
assert.equal(T.frame(model, -1).length, 0);
assert.equal(T.frame(model, 0).length, 1);
assert.deepEqual(
  T.frame(model, 5).map((s) => s.track.id),
  [1, 2],
  'Parent must disappear exactly when children are born',
);
assert.deepEqual(
  T.frame(model, 7).map((s) => s.track.id),
  [1],
  'Ablated particles must not remain on the map',
);
const final = T.frame(model, 15);
assert.equal(final.length, 1);
assert.equal(final[0].landed, true);
assert.equal(final[0].row[2], 5);
assert.ok(
  Math.abs(T.frame(model, 10)[0].row[6] - 180) < 1e-9,
  'Interpolation must cross the dateline via the shortest path',
);
assert.ok(model.bounds[1][1] - model.bounds[0][1] < 3);
assert.ok(model.markers.some((m) => m.kind === 'height' && m.row[1] === 60));
assert.ok(model.markers.some((m) => m.kind === 'dark' && Math.abs(m.row[7] - 3) < 1e-9));
assert.equal(model.markers.find((m) => m.kind === 'breakup').pressure, 4e6);
assert.throws(() => T.prepare(flight, values, 5), /coordinates unavailable/);
assert.throws(() => T.prepare({ verified: false }, values, 8), /unavailable/);
for (let id = 3; id <= 1002; id++)
  add(id, 0, 1, 'ground', rows(5, 30, 1 / id, 179.5), rows(15, 0, 0.1 / id, 179.6));
const limited = prepare(120);
assert.equal(limited.leaves, 120);
assert.equal(limited.tracks.length, 121);
assert.equal(T.frame(limited, 10).length, 120);
assert.ok(limited.tracks.some((t) => t.id === 0));
assert.ok(limited.tracks.some((t) => t.id === 1));
assert.deepEqual(
  limited.tracks.map((t) => t.id),
  prepare(120).tracks.map((t) => t.id),
  'Display sampling must be reproducible',
);
assert.equal(prepare(1).leaves, 1);
const origin = [0, 100, 1, 10, 0, 0, 0, 10],
  destination = [10, 1.3, 0.1, 0.1, 0, 0, 1, 0.1];
const profileTracks = [{ generation: 0, rows: [origin, destination] }];
const east = T.profileGeometry(profileTracks, { velocity_ecef_mps: [-10, 1, 0] });
assert.ok(east.hasHeading);
assert.equal(east.east, 1);
assert.equal(east.north, 0);
assert.ok(Math.abs(T.projectProfile(east, destination).x - 111.19508) < 1e-4);
assert.equal(
  T.projectProfile(east, destination).y,
  1.3,
  'Ground altitude must remain WGS84, not be clamped to zero',
);
assert.ok(
  T.projectProfile(east, [...destination.slice(0, 6), -1, 0.1]).x < 0,
  'Backward wind drift must retain its sign',
);
assert.equal(T.projectProfile(east, [0, 0, 1, 0, 0, 1, 0, 0]).x, 0, 'Lateral displacement is projected out');
const north = T.profileGeometry(profileTracks, { velocity_ecef_mps: [-10, 0, 1] });
assert.ok(Math.abs(T.projectProfile(north, [0, 0, 1, 0, 0, 1, 0, 0]).x - 111.19508) < 1e-4);
const vertical = T.profileGeometry([{ generation: 0, rows: [origin] }], { velocity_ecef_mps: [-10, 0, 0] });
assert.equal(
  vertical.hasHeading,
  false,
  'Vertical entry has no heading; label the northward fallback explicitly',
);
assert.ok(vertical.xmax > vertical.xmin && vertical.ymax > vertical.ymin);
assert.ok(
  T.projectProfile(model.profile, model.tracks[1].rows.at(-1)).x > 0,
  'Dateline crossing must preserve signed along-track position',
);
const context = new Proxy(
  {},
  {
    get: (target, key) => target[key] ?? (() => {}),
    set: (target, key, value) => ((target[key] = value), true),
  },
);
for (const time of [0, 5, 7, 15]) {
  const particles = T.frame(model, time),
    view = T.drawProfile(context, model, time, 360, 400, () => [10, 50, 100], particles);
  assert.equal(
    view.particles,
    particles,
    'The side projection must draw the exact shared frame, not resample it',
  );
  for (const state of particles) assert.ok(Object.values(view.point(state.row)).every(Number.isFinite));
}
if (process.argv[2]) {
  const fs = require('node:fs'),
    html = fs.readFileSync(process.argv[2], 'utf8'),
    R = require('../output/web/report.js');
  const data = JSON.parse(
    html.match(/<script id="field-data" type="application\/json">([\s\S]*?)<\/script>/)[1],
  );
  for (const run of data.runs) {
    const marker = `<script id="${run.flight.element}" type="application/json">`,
      start = html.indexOf(marker) + marker.length;
    const f = JSON.parse(html.slice(start, html.indexOf('</script>', start)));
    const v = R.floats(
        require('node:zlib').inflateSync(Buffer.from(f.samples, 'base64')).toString('base64'),
        true,
      ),
      m = T.prepare(f, v, data.trace_columns.length, 120, 3, run.entry);
    for (const particle of T.frame(m, m.end)) {
      const p = data.impacts.find((p) => p[5] === run.run_id && p[6] === particle.track.id);
      assert.ok(p);
      assert.ok(Math.abs(p[0] - particle.row[5]) < 1e-7);
      assert.ok(Math.abs(p[1] - particle.row[6]) < 1e-7);
      assert.equal(p[2], particle.row[2]);
      assert.equal(T.projectProfile(m.profile, particle.row).y, particle.track.last[1]);
    }
    for (let i = 0; i < 100; i++) assert.ok(T.frame(m, m.start + ((m.end - m.start) * i) / 99).length <= 120);
  }
  console.log(`Animation ground endpoints verified against all ${data.runs.length} archived runs.`);
}
console.log(
  'Trajectory birth/death, cutoff, dateline, bounded sampling and synchronized signed profile tests passed.',
);
