(function (root) {
  'use strict';
  const deltaLongitude = (a, b) => ((((b - a + 540) % 360) + 360) % 360) - 180;
  function interpolate(a, b, f) {
    const row = a.map((v, i) => v + (b[i] - v) * f);
    row[6] = a[6] + deltaLongitude(a[6], b[6]) * f;
    return row;
  }
  function at(track, time) {
    const rows = track.rows,
      last = rows.at(-1);
    if (time < rows[0][0]) return null;
    if (time >= last[0]) return track.end_event === 'ground' ? { row: last, landed: true } : null;
    let low = 0,
      high = rows.length - 1;
    while (high - low > 1) {
      const mid = (low + high) >> 1;
      if (rows[mid][0] <= time) low = mid;
      else high = mid;
    }
    const a = rows[low],
      b = rows[high];
    return { row: interpolate(a, b, (time - a[0]) / (b[0] - a[0])), landed: false };
  }
  function crossing(rows, column, value) {
    for (let i = 1; i < rows.length; i++)
      if (rows[i - 1][column] >= value && rows[i][column] < value)
        return interpolate(
          rows[i - 1],
          rows[i],
          (rows[i - 1][column] - value) / (rows[i - 1][column] - rows[i][column]),
        );
    return null;
  }
  function profileGeometry(tracks, entry) {
    const first = (tracks.find((t) => t.generation === 0) ?? tracks[0]).rows[0],
      rad = Math.PI / 180;
    const latitude = first[5] * rad,
      longitude = first[6] * rad,
      sinLat = Math.sin(latitude),
      cosLat = Math.cos(latitude);
    const v = entry?.velocity_ecef_mps;
    let east = 0,
      north = 0;
    if (v?.length === 3 && v.every(Number.isFinite)) {
      east = -Math.sin(longitude) * v[0] + Math.cos(longitude) * v[1];
      north = -sinLat * Math.cos(longitude) * v[0] - sinLat * Math.sin(longitude) * v[1] + cosLat * v[2];
    } else {
      for (const row of (tracks.find((t) => t.generation === 0) ?? tracks[0]).rows) {
        const lat = row[5] * rad,
          dl = deltaLongitude(first[6], row[6]) * rad;
        east = Math.cos(lat) * Math.sin(dl);
        north = cosLat * Math.sin(lat) - sinLat * Math.cos(lat) * Math.cos(dl);
        if (Math.hypot(east, north) > 1e-10) break;
      }
    }
    const norm = Math.hypot(east, north),
      hasHeading = norm > 1e-10;
    east = hasHeading ? east / norm : 0;
    north = hasHeading ? north / norm : 1;
    const geometry = { latitude: first[5], longitude: first[6], east, north, hasHeading };
    let xmin = 0,
      xmax = 0,
      ymin = 0,
      ymax = 0;
    for (const track of tracks)
      for (const row of track.rows) {
        const p = projectProfile(geometry, row);
        xmin = Math.min(xmin, p.x);
        xmax = Math.max(xmax, p.x);
        ymin = Math.min(ymin, p.y);
        ymax = Math.max(ymax, p.y);
      }
    const dx = Math.max(1, xmax - xmin),
      dy = Math.max(1, ymax - ymin);
    return {
      ...geometry,
      xmin: xmin - 0.04 * dx,
      xmax: xmax + 0.04 * dx,
      ymin: ymin - 0.04 * dy,
      ymax: ymax + 0.04 * dy,
    };
  }
  function projectProfile(geometry, row) {
    const rad = Math.PI / 180,
      lat = row[5] * rad,
      lat0 = geometry.latitude * rad,
      dl = deltaLongitude(geometry.longitude, row[6]) * rad;
    const east = Math.cos(lat) * Math.sin(dl),
      north = Math.cos(lat0) * Math.sin(lat) - Math.sin(lat0) * Math.cos(lat) * Math.cos(dl);
    const dot = Math.sin(lat0) * Math.sin(lat) + Math.cos(lat0) * Math.cos(lat) * Math.cos(dl);
    // Signed great-circle along-track distance; lateral wind drift is projected out.
    return { x: 6371.0088 * Math.atan2(east * geometry.east + north * geometry.north, dot), y: row[1] };
  }
  function prepare(flight, values, stride, limit = 120, threshold = 3, entry = null) {
    if (!flight?.verified || stride < 7)
      throw new Error('Trajectory coordinates unavailable. Regenerate this report from the saved run.');
    limit = Math.max(1, Math.floor(limit));
    const byId = new Map(flight.traces.map((t) => [t.id, t])),
      parents = new Set(flight.traces.filter((t) => t.parent_id !== t.id).map((t) => t.parent_id));
    const leaves = flight.traces
      .filter((t) => !parents.has(t.id))
      .sort((a, b) => b.last[2] - a.last[2] || a.id - b.id);
    const chosen = new Set(),
      count = Math.min(limit, leaves.length);
    for (let i = 0; i < count; i++) {
      let t = leaves[count === 1 ? 0 : Math.round((i * (leaves.length - 1)) / (count - 1))];
      while (t && !chosen.has(t.id)) {
        chosen.add(t.id);
        t = byId.get(t.parent_id);
      }
    }
    const tracks = flight.traces
      .filter((t) => chosen.has(t.id))
      .map((t) => {
        const rows = [];
        for (let i = 0; i < t.count; i++)
          rows.push(Array.from(values.subarray(t.offset + i * stride, t.offset + (i + 1) * stride)));
        rows[0] = [...t.first];
        rows[rows.length - 1] = [...t.last];
        if (
          !rows.length ||
          rows.some((r, i) => r.length < 7 || !r.every(Number.isFinite) || (i > 0 && r[0] < rows[i - 1][0]))
        )
          throw new Error('Invalid trajectory samples');
        for (let i = 1; i < rows.length; i++)
          rows[i][6] = rows[i - 1][6] + deltaLongitude(rows[i - 1][6], rows[i][6]);
        return { ...t, rows };
      });
    const children = new Map();
    for (const t of tracks)
      if (t.id !== t.parent_id) {
        if (!children.has(t.parent_id)) children.set(t.parent_id, []);
        children.get(t.parent_id).push(t);
      }
    const main = [],
      visited = new Set();
    let branch = tracks.find((t) => t.generation === 0);
    while (branch && !visited.has(branch.id)) {
      visited.add(branch.id);
      main.push(...branch.rows);
      branch = children
        .get(branch.id)
        ?.slice()
        .sort((a, b) => b.first[2] - a.first[2] || a.id - b.id)[0];
    }
    const markers = [];
    for (let height = 10; height <= Math.min(150, main[0]?.[1] ?? 0); height += 10) {
      const row = crossing(main, 1, height);
      if (row) markers.push({ kind: 'height', row, label: height + ' km' });
    }
    for (const b of flight.breakups ?? [])
      if (chosen.has(b.fragment_id) && b.state.length >= 7)
        markers.push({
          kind: 'breakup',
          row: [...b.state],
          label: 'Breakup #' + b.fragment_id,
          pressure: b.pressure_pa,
          fragmentId: b.fragment_id,
        });
    const dark = stride >= 8 ? crossing(main, 7, threshold) : null;
    if (dark) markers.push({ kind: 'dark', row: dark, label: 'Dark-flight threshold (displayed branch)' });
    let start = Infinity,
      end = 0,
      south = 90,
      north = -90,
      west = Infinity,
      east = -Infinity,
      lo = Infinity,
      hi = -Infinity;
    const origin = tracks[0]?.rows[0][6] ?? 0;
    for (const t of tracks)
      for (const row of t.rows) {
        start = Math.min(start, row[0]);
        end = Math.max(end, row[0]);
        south = Math.min(south, row[5]);
        north = Math.max(north, row[5]);
        row[6] = origin + deltaLongitude(origin, row[6]);
        west = Math.min(west, row[6]);
        east = Math.max(east, row[6]);
        if (row[2] > 0) {
          lo = Math.min(lo, Math.log10(row[2]));
          hi = Math.max(hi, Math.log10(row[2]));
        }
      }
    if (!tracks.length) throw new Error('No trajectory samples available.');
    for (const marker of markers) marker.row[6] = origin + deltaLongitude(origin, marker.row[6]);
    for (const track of flight.traces) {
      start = Math.min(start, track.first[0]);
      end = Math.max(end, track.last[0]);
    }
    return {
      tracks,
      markers,
      main,
      start,
      end,
      lo,
      hi,
      profile: profileGeometry(tracks, entry),
      bounds: [
        [south, west],
        [north, east],
      ],
      total: flight.trace_selection?.total ?? flight.traces.length,
      leaves: count,
      limit,
      firstBreakup: (flight.breakups ?? []).reduce((first, b) => Math.min(first, b.state[0]), end),
    };
  }
  function frame(model, time) {
    return model.tracks
      .map((track) => {
        const state = at(track, time);
        return state ? { ...state, track } : null;
      })
      .filter(Boolean);
  }
  function drawScene(ctx, model, time, point, palette, particles, light = false) {
    ctx.save();
    for (const track of model.tracks) {
      ctx.beginPath();
      track.rows.forEach((r, i) => {
        const p = point(r);
        if (i) ctx.lineTo(p.x, p.y);
        else ctx.moveTo(p.x, p.y);
      });
      ctx.strokeStyle =
        track.generation === 0
          ? light
            ? '#85652f'
            : '#fff2bc'
          : light
            ? 'rgba(76,103,116,.18)'
            : 'rgba(231,247,252,.32)';
      ctx.lineWidth = track.generation === 0 ? 3 : 1;
      ctx.stroke();
    }
    for (const state of particles) {
      const row = state.row,
        p = point(row),
        fraction = (Math.log10(Math.max(1e-30, row[2])) - model.lo) / Math.max(1e-9, model.hi - model.lo);
      const radius = 3 + 3 * Math.max(0, Math.min(1, fraction));
      if (!state.landed) {
        const previous = at(state.track, Math.max(state.track.rows[0][0], time - 1.5));
        if (previous) {
          const a = point(previous.row);
          ctx.beginPath();
          ctx.moveTo(a.x, a.y);
          ctx.lineTo(p.x, p.y);
          ctx.strokeStyle = light ? '#496878' : '#fff9';
          ctx.lineWidth = 2;
          ctx.stroke();
        }
      }
      ctx.beginPath();
      ctx.arc(p.x, p.y, radius, 0, 2 * Math.PI);
      ctx.fillStyle = `rgb(${palette(fraction).join(',')})`;
      ctx.fill();
      ctx.strokeStyle = state.landed || light ? '#213f48' : '#fff';
      ctx.lineWidth = state.landed ? 1 : 1.5;
      ctx.stroke();
    }
    ctx.restore();
    return particles;
  }
  function draw(ctx, model, time, point, palette) {
    return drawScene(ctx, model, time, (r) => point(r[5], r[6]), palette, frame(model, time));
  }
  function profileTicks(min, max, count) {
    const raw = (max - min) / count,
      base = 10 ** Math.floor(Math.log10(raw)),
      scaled = raw / base;
    const step = (scaled <= 1 ? 1 : scaled <= 2 ? 2 : scaled <= 5 ? 5 : 10) * base,
      values = [];
    for (let i = Math.ceil(min / step); i * step <= max; i++) values.push(i * step);
    return { values, digits: Math.max(0, -Math.floor(Math.log10(step))) };
  }
  function drawProfile(ctx, model, time, width, height, palette, particles = frame(model, time)) {
    const g = model.profile,
      left = 56,
      right = width - 18,
      top = 40,
      bottom = height - 54;
    const point = (row) => {
      const p = projectProfile(g, row);
      return {
        x: left + ((p.x - g.xmin) / (g.xmax - g.xmin)) * (right - left),
        y: bottom - ((p.y - g.ymin) / (g.ymax - g.ymin)) * (bottom - top),
      };
    };
    ctx.save();
    ctx.fillStyle = '#fff';
    ctx.fillRect(0, 0, width, height);
    ctx.font = '12px Arial';
    ctx.lineWidth = 1;
    ctx.textAlign = 'left';
    ctx.fillStyle = '#425d69';
    ctx.fillText('WGS84 height / km', left, 20);
    const xs = profileTicks(g.xmin, g.xmax, width < 450 ? 4 : 6),
      ys = profileTicks(g.ymin, g.ymax, 6);
    for (const v of xs.values) {
      const x = left + ((v - g.xmin) / (g.xmax - g.xmin)) * (right - left);
      ctx.beginPath();
      ctx.moveTo(x, top);
      ctx.lineTo(x, bottom);
      ctx.strokeStyle = '#e0e7eb';
      ctx.stroke();
      ctx.textAlign = 'center';
      ctx.fillStyle = '#425d69';
      ctx.fillText(v.toFixed(xs.digits), x, bottom + 18);
    }
    for (const v of ys.values) {
      const y = bottom - ((v - g.ymin) / (g.ymax - g.ymin)) * (bottom - top);
      ctx.beginPath();
      ctx.moveTo(left, y);
      ctx.lineTo(right, y);
      ctx.strokeStyle = '#e0e7eb';
      ctx.stroke();
      ctx.textAlign = 'right';
      ctx.fillStyle = '#425d69';
      ctx.fillText(v.toFixed(ys.digits), left - 9, y + 4);
    }
    ctx.strokeStyle = '#80959f';
    ctx.beginPath();
    ctx.moveTo(left, top);
    ctx.lineTo(left, bottom);
    ctx.lineTo(right, bottom);
    ctx.stroke();
    ctx.textAlign = 'center';
    ctx.fillStyle = '#425d69';
    ctx.fillText(
      g.hasHeading ? 'Along-entry distance / km' : 'Northward distance / km',
      (left + right) / 2,
      height - 13,
    );
    ctx.save();
    ctx.beginPath();
    ctx.rect(left, top, right - left, bottom - top);
    ctx.clip();
    drawScene(ctx, model, time, point, palette, particles, true);
    for (const marker of model.markers)
      if (marker.kind !== 'height') {
        const p = point(marker.row);
        ctx.strokeStyle =
          marker.row[0] <= time ? (marker.kind === 'breakup' ? '#c95e23' : '#00808b') : '#9bacb3';
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        if (marker.kind === 'breakup') {
          ctx.moveTo(p.x - 4, p.y - 4);
          ctx.lineTo(p.x + 4, p.y + 4);
          ctx.moveTo(p.x - 4, p.y + 4);
          ctx.lineTo(p.x + 4, p.y - 4);
        } else ctx.arc(p.x, p.y, 4, 0, 2 * Math.PI);
        ctx.stroke();
      }
    ctx.restore();
    ctx.restore();
    return { point, particles };
  }
  const api = { prepare, at, frame, crossing, draw, profileGeometry, projectProfile, drawProfile };
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.StrewnTrajectory = api;
})(globalThis);
