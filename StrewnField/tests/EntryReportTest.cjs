const fs = require('node:fs'),
  assert = require('node:assert/strict'),
  path = require('node:path');
const file = process.argv[2],
  html = fs.readFileSync(file, 'utf8');
const data = JSON.parse(
  html.match(/<script id="field-data" type="application\/json">([\s\S]*?)<\/script>/)[1],
);
const input = process.argv[3] ?? path.resolve(file, '../..');
const entries = JSON.parse(fs.readFileSync(path.join(input, 'entries.json'), 'utf8')).entries;
const summary = data.entry_summary;
assert.ok(summary);
const weighted = entries.filter((e) => e.weight > 0),
  weight = weighted.reduce((s, e) => s + e.weight, 0);
assert.equal(summary.total_entries, entries.length);
assert.equal(summary.positive_weight_entries, weighted.length);
const epoch = weighted[0].epoch_mjd_utc,
  dt = weighted.reduce((s, e) => s + e.weight * (e.epoch_mjd_utc - epoch) * 86400, 0) / weight;
assert.ok(Math.abs(summary.epoch_mjd_utc - epoch - dt / 86400) < 1e-11);
const speeds = weighted.map((e) => Math.hypot(...e.velocity_ecef_mps) / 1000),
  speedMean = speeds.reduce((s, v, i) => s + v * weighted[i].weight, 0) / weight;
assert.ok(Math.abs(summary.ecef_speed_km_s.mean - speedMean) < 1e-11);
if (weighted.length > 1) {
  const sigma = Math.sqrt(
    weighted.reduce((s, e) => s + e.weight * ((e.epoch_mjd_utc - epoch) * 86400 - dt) ** 2, 0) / weight,
  );
  assert.ok(Math.abs(summary.time_sigma_s - sigma) < 1e-10);
  const speedSigma = Math.sqrt(
    speeds.reduce((s, v, i) => s + weighted[i].weight * (v - speedMean) ** 2, 0) / weight,
  );
  assert.ok(Math.abs(summary.ecef_speed_km_s.sigma - speedSigma) < 1e-11);
} else {
  assert.equal(summary.time_sigma_s, null);
  assert.equal(summary.ecef_speed_km_s.sigma, null);
  assert.equal(summary.ellipse.semi_major_m, null);
}
assert.ok(html.includes('id="entry-statistics"'));
assert.ok(html.includes("page('Atmospheric entry'"));
console.log(
  `Original-entry statistics verified: ${weighted.length} entries, ${data.runs.length} atmospheric realizations.`,
);
