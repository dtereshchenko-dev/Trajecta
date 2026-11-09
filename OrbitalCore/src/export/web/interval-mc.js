(function (root) {
  'use strict';
  const numeric = (v) => typeof v === 'number' && Number.isFinite(v);
  function utc(mjd) {
    if (!numeric(mjd)) return 'UTC unavailable';
    const day = Math.floor(mjd),
      ms = Math.round((mjd - day) * 86400000);
    return new Date((day - 40587) * 86400000 + ms).toISOString().replace('T', ' ').replace('Z', '');
  }
  function longitudeFrame(points) {
    const values = points.map((p) => ((p.longitude_deg % 360) + 360) % 360).sort((a, b) => a - b);
    if (!values.length) return 0;
    let gap = -1,
      start = values[0];
    values.forEach((value, i) => {
      const next = values[(i + 1) % values.length] + (i === values.length - 1 ? 360 : 0);
      if (next - value > gap) {
        gap = next - value;
        start = next % 360;
      }
    });
    const centre = start + (360 - gap) / 2;
    return centre > 180 ? centre - 360 : centre;
  }
  const unwrap = (longitude, centre) => centre + ((((longitude - centre + 180) % 360) + 360) % 360) - 180;
  const ordered = (entries, event, order) =>
    entries
      .filter((p) => event === 'all' || p.event === Number(event))
      .slice()
      .sort((a, b) =>
        order === 'clone' ? a.clone - b.clone || a.epoch - b.epoch : a.epoch - b.epoch || a.clone - b.clone,
      );
  function mapProfile(target, points) {
    if (target === 'Moon' || target === 'Mars') {
      const moon = target === 'Moon';
      return {
        kind: 'wmts',
        geographic: true,
        radius: moon ? 1737400 : 3396190,
        maxNativeZoom: moon ? 8 : 7,
        name: moon ? 'LRO WAC' : 'Viking MDIM 2.1',
        url: `https://trek.nasa.gov/tiles/${target}/EQ/${
          moon ? 'LRO_WAC_Mosaic_Global_303ppd_v02' : 'Mars_Viking_MDIM21_ClrMosaic_global_232m'
        }/1.0.0/default/default028mm/{z}/{y}/{x}.jpg`,
        attribution: moon
          ? 'LRO WAC &copy; NASA / GSFC / ASU / USGS'
          : 'Viking MDIM 2.1 &copy; NASA / JPL / USGS / NASA Ames',
      };
    }
    if (target === 'Earth' && !points.some((p) => Math.abs(p.latitude_deg) > 85))
      return {
        kind: 'xyz',
        geographic: false,
        radius: 6371000,
        maxNativeZoom: 19,
        name: 'Esri World Imagery',
        url: 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}?blankTile=false',
        attribution: 'Imagery &copy; Esri, Vantor, Earthstar Geographics, GIS User Community',
      };
    return { kind: 'grid', geographic: true, name: 'Geographic grid' };
  }
  function imageryLayer(L, profile, overview = false) {
    const options = { maxZoom: 22, maxNativeZoom: profile.maxNativeZoom, attribution: profile.attribution };
    if (profile.kind === 'xyz') return L.tileLayer(profile.url, options);
    if (profile.kind === 'wmts')
      return L.tileLayer(profile.url, {
        ...options,
        minNativeZoom: 0,
        maxNativeZoom: overview ? 0 : profile.maxNativeZoom,
        bounds: [
          [-90, -180],
          [90, 180],
        ],
        zIndex: overview ? 1 : 2,
        className: overview ? 'planet-overview' : 'planet-detail',
        updateWhenIdle: true,
        updateWhenZooming: false,
        keepBuffer: 4,
      });
    return null;
  }
  const api = { utc, longitudeFrame, unwrap, ordered, mapProfile, imageryLayer };
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = api;
    return;
  }
  root.OrbitalReportMap = api;
  root.addEventListener('beforeprint', () => {
    const frame = document.getElementById('map-print-frame'),
      map = document.getElementById('map');
    if (!frame || !map || !map.offsetWidth || !map.offsetHeight) return;
    const width = map.offsetWidth,
      height = map.offsetHeight;
    const scale = Math.min(1, (260 * 96) / 25.4 / width, (105 * 96) / 25.4 / height);
    frame.style.setProperty('--print-map-width', `${width}px`);
    frame.style.setProperty('--print-map-height', `${height}px`);
    frame.style.setProperty('--print-map-scale', String(scale));
    frame.style.setProperty('--print-frame-width', `${width * scale}px`);
    frame.style.setProperty('--print-frame-height', `${height * scale}px`);
  });
  if (!document.getElementById('mc-data')) return;
  const el = (id) => document.getElementById(id),
    data = JSON.parse(el('mc-data').textContent);
  const fmt = (v, digits = 7) =>
    numeric(v) ? Number(v).toLocaleString('en-US', { maximumSignificantDigits: digits }) : '\u2014';
  const fixed = (v, digits = 6) => (numeric(v) ? v.toFixed(digits) : '\u2014');
  const target = data.target_name;
  const earth = target === 'Earth',
    height = data.altitude_km;
  const coordinateModel = data.body_fixed_frame
    ? `${data.body_fixed_frame}; ${data.body_fixed_model}`
    : earth
      ? data.earth_frame || 'Earth frame not recorded'
      : 'Legacy approximate orientation; frame not recorded';
  const boundary =
    earth && height > 0
      ? `Atmospheric entry at ${fmt(height)} km`
      : height === 0
        ? `${target} reference-surface crossings`
        : `${target} boundary crossings at ${fmt(height)} km`;
  const points = data.entries,
    events = data.events,
    pageSize = 100;
  const profile = mapProfile(target, points);
  let page = 0,
    selected = null,
    map = null,
    layer = null,
    imagery = null,
    overview = null,
    centre = 0,
    markers = new Map();
  function table(id, caption, headers, rows, onRow) {
    const t = document.createElement('table');
    t.createCaption().textContent = caption;
    const head = t.createTHead().insertRow();
    headers.forEach((title) => {
      const th = document.createElement('th');
      th.scope = 'col';
      th.textContent = title;
      head.append(th);
    });
    const body = t.createTBody();
    rows.forEach((values, i) => {
      const row = body.insertRow();
      values.forEach((value) => {
        row.insertCell().textContent = value;
      });
      if (onRow) onRow(row, i);
    });
    el(id).replaceChildren(t);
  }
  el('subtitle').textContent =
    `${data.method === 'covariance' ? 'Covariance MC' : 'Observation MC'} | ${target} | ${boundary}`;
  const metrics = [
    ['Completed trials', `${data.completed} / ${data.requested}`],
    ['Clones with a crossing', fmt(data.entry_clones)],
    ['Crossing events', events.length],
    ['Failed trials', fmt(data.failed_fit + data.failed_propagation)],
  ];
  metrics.forEach(([label, value]) => {
    const m = document.createElement('div');
    m.className = 'metric';
    const caption = document.createElement('div');
    caption.textContent = label;
    const v = document.createElement('div');
    v.className = 'value';
    v.textContent = value;
    m.append(caption, v);
    el('metrics').append(m);
  });
  const warnings = [];
  if (data.cancelled)
    warnings.push('Cancelled run: completed trials are retained; partial trajectories are excluded.');
  if (data.completed !== data.requested)
    warnings.push(
      `Incomplete sample: ${data.failed_fit} failed fits, ${data.failed_propagation} failed propagations, ${data.cancelled_clones} cancelled, ${data.unstarted} not started. Missing trials can bias probabilities.`,
    );
  if (data.nominal.status !== 1) warnings.push(`Nominal propagation unavailable: ${data.nominal.error}`);
  if (data.entries_before_interval)
    warnings.push(
      `${data.entries_before_interval} clones crossed the boundary before this interval; no later crossings are inferred.`,
    );
  el('warnings').textContent = warnings.join('\n');
  el('statistics-note').textContent =
    `${utc(data.start_utc)} to ${utc(data.end_utc)} UTC. Seed ${data.seed}; ${data.integrator}; ${data.workers} workers. ` +
    'Probabilities use all complete trials, including those without a crossing. Intervals are Wilson 95%, conditional on successful propagation. Zero entries does not establish zero risk.';
  table(
    'event-table',
    'Physical boundary crossings',
    ['Event', 'First crossing, UTC', 'Last crossing, UTC', 'Entries', 'P', '95% lower', '95% upper'],
    events.map((e, i) => [
      i + 1,
      utc(e.start_utc),
      utc(e.end_utc),
      e.entries,
      fmt(e.probability),
      fmt(e.probability_low),
      fmt(e.probability_high),
    ]),
  );
  el('event-note').textContent =
    `Only physical boundary crossings are listed. Gaps over ${fmt(data.event_max_gap_days * 24)} hours separate groups. ` +
    'Grouping is for presentation only, not a search window or proof of distinct dynamical returns. Nominal is not counted as an MC trial. ' +
    (events.length ? '' : 'No boundary crossings in completed trials. This does not establish zero risk.');
  const config = { ...data };
  delete config.entries;
  delete config.events;
  delete config.nominal;
  el('configuration').textContent = JSON.stringify(config, null, 2);
  ['map-title', 'entry-title'].forEach((id) => {
    el(id).textContent = boundary;
  });
  ['map-event', 'entry-filter'].forEach((id) => {
    el(id).add(new Option('All events', 'all'));
    events.forEach((e, i) =>
      el(id).add(new Option(`${i + 1}: ${utc(e.start_utc)} UTC (${e.entries} crossings)`, String(i))),
    );
  });
  el('entry-note').textContent =
    `Coordinate frame: ${coordinateModel}. ` +
    'Latitude is positive north, longitude positive east. Earth coordinates are geodetic on WGS84; other targets use their reference surface. ' +
    'Speed is target-relative inertial speed, not atmosphere-relative speed. Nominal is separate from MC trials. ' +
    (height > 0
      ? 'These points are crossings of the selected altitude, not predicted ground-fall locations. Atmospheric descent is not modelled here.'
      : 'The reference surface is not local terrain; atmospheric drag and fragmentation are not modelled here.');
  function rows() {
    return ordered(points, el('entry-filter').value, el('entry-order').value);
  }
  function renderTable() {
    const all = rows(),
      pages = Math.max(1, Math.ceil(all.length / pageSize));
    page = Math.max(0, Math.min(page, pages - 1));
    const slice = all.slice(page * pageSize, (page + 1) * pageSize);
    table(
      'entry-table',
      boundary,
      ['Clone', 'Event', 'Date / time UTC', 'Latitude, deg N', 'Longitude, deg E', 'Height, km', 'km/s', ''],
      slice.map((p) => [
        p.clone < 0 ? 'Nominal' : p.clone + 1,
        p.event === null ? '\u2014' : p.event + 1,
        utc(p.epoch_utc),
        fixed(p.latitude_deg),
        fixed(p.longitude_deg),
        fmt(height),
        fmt(Math.hypot(...p.relative_velocity_km_s)),
        '',
      ]),
      (row, i) => {
        const p = slice[i];
        if (p === selected) row.className = 'selected';
        const button = document.createElement('button');
        button.textContent = 'Map';
        button.title = 'Show crossing on map';
        button.onclick = () => focusPoint(p);
        row.cells[7].append(button);
      },
    );
    el('entry-count').textContent =
      `${all.filter((p) => p.clone >= 0).length} MC crossings; ${all.filter((p) => p.clone < 0).length} nominal`;
    el('page-count').textContent = `${page + 1} / ${pages}`;
    el('previous').disabled = page === 0;
    el('next').disabled = page + 1 === pages;
  }
  function popup(p) {
    const node = document.createElement('div');
    node.textContent =
      `${p.clone < 0 ? 'Nominal' : `Clone ${p.clone + 1}`} | ${utc(p.epoch_utc)} UTC | ` +
      `${fixed(p.latitude_deg)} deg N, ${fixed(p.longitude_deg)} deg E | h = ${fmt(height)} km`;
    return node;
  }
  function visiblePoints() {
    return ordered(points, el('map-event').value, 'time').filter(
      (p) => el('show-nominal').checked || p.clone >= 0,
    );
  }
  function fitMap() {
    const bounds = visiblePoints().map((p) => [p.latitude_deg, unwrap(p.longitude_deg, centre)]);
    if (bounds.length) map.fitBounds(bounds, { padding: [40, 40], maxZoom: 13, animate: false });
    else map.setView([0, 0], 2);
  }
  function renderMap(fit = true) {
    if (!map) return;
    layer.clearLayers();
    markers.clear();
    const visible = visiblePoints();
    centre = longitudeFrame(visiblePoints());
    visible.forEach((p) => {
      const marker = L.circleMarker([p.latitude_deg, unwrap(p.longitude_deg, centre)], {
        radius: p.clone < 0 ? 6 : 3.5,
        color: p.clone < 0 ? '#4f3515' : '#702925',
        fillColor: p.clone < 0 ? '#ffca68' : '#ff5b57',
        weight: p.clone < 0 ? 2 : 0.7,
        fillOpacity: 0.9,
      })
        .bindPopup(popup(p))
        .addTo(layer);
      marker.on('click', () => {
        selected = p;
        renderTable();
      });
      markers.set(p, marker);
    });
    el('map-empty').textContent = visiblePoints().length
      ? ''
      : 'No boundary crossings in the selected event.';
    if (fit) fitMap();
  }
  function initMap() {
    if (map) return;
    const crs = profile.geographic
      ? L.extend({}, L.CRS.EPSG4326, {
          distance: (a, b) => (L.CRS.Earth.distance(a, b) * (profile.radius || 6371000)) / 6371000,
        })
      : L.CRS.EPSG3857;
    map = L.map('map', {
      preferCanvas: true,
      zoomAnimation: false,
      crs,
      minZoom: 0,
      maxZoom: 22,
      wheelDebounceTime: 120,
      wheelPxPerZoomLevel: 100,
    });
    map.on('resize', () => map.closePopup());
    map.setView([0, 0], 2);
    layer = L.layerGroup().addTo(map);
    if (profile.radius) L.control.scale({ imperial: false }).addTo(map);
    imagery = imageryLayer(L, profile);
    overview = profile.kind === 'wmts' ? imageryLayer(L, profile, true) : null;
    el('imagery').disabled = !imagery;
    if (!imagery) el('imagery').checked = false;
    el('imagery-label').textContent = earth ? 'Satellite imagery' : 'Surface imagery';
    el('imagery').title = profile.name;
    const note =
      'Coordinates and markers are embedded and work offline. ' +
      (imagery
        ? `${profile.name}${overview ? ' via NASA Trek' : ''}; optional online imagery. `
        : 'No imagery available for this target or polar extent. ') +
      (profile.kind === 'wmts' ? 'Planetocentric latitude, east-positive longitude. ' : '') +
      (!earth && !data.body_fixed_frame
        ? 'Legacy archive: approximate orientation; coordinates have not been recalculated.'
        : '');
    let tileErrors = 0;
    el('map-meta').textContent = note;
    if (imagery) {
      imagery.on('loading', () => {
        tileErrors = 0;
        el('map-meta').textContent = note;
      });
      imagery.on('tileerror', () => {
        ++tileErrors;
        el('map-meta').textContent =
          note + 'Some imagery is unavailable; coordinates and markers remain visible.';
      });
      imagery.on('load', () => {
        if (!tileErrors) el('map-meta').textContent = note;
      });
      el('imagery').onchange = () => {
        el('map-meta').textContent = note;
        if (el('imagery').checked) {
          overview?.addTo(map);
          imagery.addTo(map);
        } else {
          imagery.remove();
          overview?.remove();
        }
      };
      if (el('imagery').checked) {
        overview?.addTo(map);
        imagery.addTo(map);
      }
    }
    const pane = map.createPane('coordinates');
    pane.style.zIndex = 450;
    pane.style.pointerEvents = 'none';
    const canvas = L.DomUtil.create('canvas', 'coordinate-canvas');
    pane.append(canvas);
    function grid() {
      canvas.style.display = el('coordinates').checked ? 'block' : 'none';
      if (!el('coordinates').checked) return;
      const size = map.getSize(),
        dpr = window.devicePixelRatio || 1;
      canvas.width = size.x * dpr;
      canvas.height = size.y * dpr;
      canvas.style.width = size.x + 'px';
      canvas.style.height = size.y + 'px';
      L.DomUtil.setPosition(canvas, map.containerPointToLayerPoint([0, 0]));
      const ctx = canvas.getContext('2d');
      ctx.scale(dpr, dpr);
      StrewnGraticule.draw(
        ctx,
        size.x,
        size.y,
        (x, y) => map.containerPointToLatLng([x, y]),
        (lat, lon) => {
          const p = map.latLngToContainerPoint([lat, lon]);
          return [p.x, p.y];
        },
      );
    }
    map.on('moveend zoomend resize', grid);
    el('coordinates').onchange = grid;
    renderMap();
    grid();
  }
  function view(name) {
    document
      .querySelectorAll('[data-view]')
      .forEach((b) => b.setAttribute('aria-selected', String(b.dataset.view === name)));
    document.querySelectorAll('.view').forEach((s) => s.classList.toggle('active', s.id === name));
    if (name === 'field') {
      initMap();
      map.invalidateSize(false);
    }
  }
  function focusPoint(p) {
    selected = p;
    el('map-event').value = p.event === null ? 'all' : String(p.event);
    el('show-nominal').checked = true;
    view('field');
    renderMap(false);
    map.setView([p.latitude_deg, unwrap(p.longitude_deg, centre)], Math.max(map.getZoom(), 12), {
      animate: false,
    });
    markers.get(p)?.openPopup();
  }
  document.querySelectorAll('[data-view]').forEach((b) => {
    b.onclick = () => view(b.dataset.view);
  });
  el('entry-filter').onchange = () => {
    page = 0;
    renderTable();
  };
  el('entry-order').onchange = () => {
    page = 0;
    renderTable();
  };
  el('previous').onclick = () => {
    --page;
    renderTable();
  };
  el('next').onclick = () => {
    ++page;
    renderTable();
  };
  el('map-event').onchange = () => renderMap();
  el('show-nominal').onchange = () => renderMap();
  el('fit').onclick = fitMap;
  el('print').onclick = () => window.print();
  renderTable();
})(typeof globalThis !== 'undefined' ? globalThis : this);
