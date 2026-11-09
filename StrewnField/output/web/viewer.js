(async function () {
  'use strict';
  const data = JSON.parse(document.getElementById('field-data').textContent);
  const R = StrewnReport,
    el = (id) => document.getElementById(id),
    runs = data.runs;
  const format = (v) =>
    v == null ? 'n/a' : Number(v).toLocaleString('en-US', { maximumSignificantDigits: 6 });
  const colors = ['#087f83', '#cb782a', '#5965a7', '#be525a', '#6a8f45'];
  const config = {
    responsive: true,
    displaylogo: false,
    modeBarButtonsToRemove: ['toImage', 'lasso2d', 'select2d'],
    scrollZoom: false,
  };
  const plots = new Set(),
    printPlots = new Map();
  let selection = 0,
    decoded = null,
    currentFlight = null,
    active = 'overview',
    flightRevision = 0,
    selectionRevision = 0;
  function layout(title, x, y, extra = {}) {
    return {
      title: { text: title, font: { size: 15 }, x: 0.08 },
      margin: { l: 75, r: 25, t: 58, b: 58 },
      font: { family: 'system-ui, Segoe UI, sans-serif', size: 12, color: '#334a54' },
      paper_bgcolor: '#fff',
      plot_bgcolor: '#fff',
      xaxis: { title: { text: x }, gridcolor: '#e3e9ec', zeroline: false, automargin: true },
      yaxis: { title: { text: y }, gridcolor: '#e3e9ec', zeroline: false, automargin: true },
      showlegend: false,
      ...extra,
    };
  }
  async function plot(id, traces, options) {
    plots.add(id);
    printPlots.set(id, structuredClone({ data: traces, layout: options }));
    await Plotly.react(el(id), traces, options, config);
  }
  function table(target, caption, head, rows) {
    const node = document.createElement('table'),
      cap = node.createCaption();
    cap.textContent = caption;
    const header = node.createTHead().insertRow();
    head.forEach((text) => {
      const th = document.createElement('th');
      th.textContent = text;
      header.append(th);
    });
    const body = node.createTBody();
    rows.forEach((row) => {
      const tr = body.insertRow();
      row.forEach((text) => {
        tr.insertCell().textContent = text;
      });
    });
    el(target).replaceChildren(node);
  }
  function bar(id, title, x, values, weights, log = false, color = colors[0], histogram = null) {
    const h = histogram ?? R.histogram(values, weights, 28, log),
      axes = layout(title, x, 'Probability per bin');
    axes.xaxis.type = log ? 'log' : 'linear';
    axes.yaxis.tickformat = '.0%';
    if (h.singleValue !== undefined) {
      const value = h.singleValue,
        center = log ? Math.log10(value) : value,
        span = log ? 0.5 : id === 'counts' ? 1 : Math.max(Math.abs(value) * 0.05, 1e-6);
      Object.assign(axes.xaxis, {
        range: [center - span, center + span],
        tickvals: [value],
        ticktext: [format(value)],
      });
      Object.assign(axes.yaxis, { range: [0, 1.18], title: { text: 'Probability mass' } });
      axes.shapes = [{ type: 'line', x0: value, x1: value, y0: 0, y1: 1, line: { color, width: 3 } }];
      return plot(
        id,
        [
          {
            type: 'scatter',
            mode: 'markers+text',
            x: [value],
            y: [1],
            marker: { color, size: 10 },
            text: [`100% at ${format(value)}`],
            textposition: 'top center',
            hovertemplate: '%{x:.6g}<br>100% of samples<extra></extra>',
          },
        ],
        axes,
      );
    }
    return plot(
      id,
      [
        {
          type: 'bar',
          x: h.x,
          y: h.y,
          width: h.width,
          marker: { color },
          customdata: h.ranges,
          hovertemplate: '%{customdata[0]:.5g} to %{customdata[1]:.5g}<br>%{y:.2%}<extra></extra>',
        },
      ],
      axes,
    );
  }
  async function overview() {
    entryTable();
    const weights = runs.map((r) => r.run_weight),
      stats = R.metrics.map((m) => R.describe(runs.map(m[2]), weights));
    const cards = [0, 1, 3, 4];
    el('metrics').replaceChildren();
    cards.forEach((i) => {
      const node = document.createElement('div');
      node.className = 'metric';
      const label = document.createElement('div');
      label.textContent = `Median ${R.metrics[i][0].toLowerCase()}`;
      const value = document.createElement('div');
      value.className = 'value';
      value.textContent = `${format(stats[i].median)} ${R.metrics[i][1]}`;
      const spread = document.createElement('div');
      spread.className = 'muted';
      spread.textContent = `P05-P95: ${format(stats[i].p05)} - ${format(stats[i].p95)}`;
      node.append(label, value, spread);
      el('metrics').append(node);
    });
    el('statistics-note').textContent =
      `${runs.length} realizations; ${runs.filter((r) => r.tracked_count === 0).length} with no survivors (included). Per-realization statistics use run weights; fragment counts include multiplicity. First-breakup statistics exclude runs without breakup. Equal-weight quantiles use linear interpolation (type 7).`;
    table(
      'statistics',
      'Ensemble statistics',
      ['Metric', 'Unit', 'N', 'Mean', 'P05', 'Median', 'P95'],
      R.metrics.map((m, i) => [
        m[0],
        m[1],
        stats[i].n,
        ...['mean', 'p05', 'median', 'p95'].map((k) => format(stats[i][k])),
      ]),
    );
    await bar(
      'landed',
      'Total landed mass per realization',
      'kg',
      runs.map((r) => r.mass_balance_kg.landed),
      weights,
    );
    await bar(
      'counts',
      'Surviving fragments per realization',
      'Count including multiplicity',
      runs.map((r) => r.ground_count),
      weights,
      false,
      colors[1],
    );
    await bar(
      'largest',
      'Largest survivor per realization',
      'kg',
      runs.map((r) => r.largest_mass_kg),
      weights,
    );
    await bar(
      'breakup-height',
      'First breakup height',
      'WGS84 height / km',
      runs.map((r) => r.first_breakup_km),
      weights,
      false,
      colors[1],
    );
    await bar(
      'mass-spectrum',
      'Surviving fragment mass spectrum',
      'kg (log)',
      data.impacts.map((p) => p[2]),
      data.impacts.map((p) => p[3] * p[4]),
      true,
      colors[0],
      data.impact_histograms?.mass,
    );
    const checks = data.summary.convergence_history ?? [],
      conv = data.job.monte_carlo?.convergence;
    const traces = ['count', 'mass'].map((q, i) => ({
      x: checks.map((c) => c.runs),
      y: checks.map((c) => c['js_' + q]),
      type: 'scatter',
      mode: 'lines+markers',
      name: q,
      line: { color: colors[i] },
      connectgaps: false,
    }));
    if (conv)
      traces.push({
        x: [0, data.summary.completed],
        y: [conv.threshold, conv.threshold],
        mode: 'lines',
        name: 'Threshold',
        line: { color: colors[3], dash: 'dash' },
      });
    await plot(
      'convergence',
      traces,
      layout('JS convergence', 'Realizations', 'Base-2 divergence', {
        showlegend: true,
        legend: { orientation: 'h', y: -0.4 },
        margin: { l: 75, r: 25, t: 58, b: 100 },
      }),
    );
    await plot(
      'mass-count',
      [
        {
          x: runs.map((r) => r.ground_count),
          y: runs.map((r) => r.mass_balance_kg.landed),
          customdata: runs.map((r, i) => [r.run_id, i]),
          mode: 'markers',
          type: 'scatter',
          marker: { size: 6, color: colors[0], opacity: 0.65 },
          hovertemplate: 'Run %{customdata[0]}<br>%{x:.5g} fragments<br>%{y:.5g} kg<extra></extra>',
        },
      ],
      layout('Realization outcomes', 'Surviving fragment count', 'Total landed mass / kg'),
    );
    el('mass-count').removeAllListeners('plotly_click');
    el('mass-count').on('plotly_click', async (e) => {
      await choose(e.points[0].customdata[1]);
      await view('flight');
    });
    const keys = ['landed', 'ablated', 'dust', 'below_cutoff', 'residual', 'unresolved'];
    const total = weights.reduce((s, w) => s + w, 0);
    await plot(
      'mass-budget',
      [
        {
          type: 'bar',
          x: keys,
          y: keys.map((k) =>
            total
              ? runs.reduce(
                  (s, r) => s + (r.run_weight * r.mass_balance_kg[k]) / r.mass_balance_kg.initial,
                  0,
                ) / total
              : 0,
          ),
          marker: { color: [...colors, '#727c80'] },
        },
      ],
      layout('Mean mass budget per realization', 'Mass destination', 'Fraction of initial mass', {
        yaxis: { tickformat: '.0%', gridcolor: '#e3e9ec', range: [0, 1] },
      }),
    );
    await bar(
      'impact-speed',
      'Impact speed distribution',
      'ECEF speed / m/s',
      data.impacts.map((p) => p[8]),
      data.impacts.map((p) => p[3] * p[4]),
      false,
      colors[1],
      data.impact_histograms?.speed,
    );
    await bar(
      'impact-time',
      'Impact time distribution',
      'Time since entry / s',
      data.impacts.map((p) => p[7]),
      data.impacts.map((p) => p[3] * p[4]),
      false,
      colors[0],
      data.impact_histograms?.time,
    );
  }
  function entryTable() {
    const s = data.entry_summary;
    if (!s) {
      el('entry-note').textContent = 'Atmospheric entry summary unavailable in this archive.';
      return;
    }
    const precise = (v, sigma, fallback = 6) =>
      v == null
        ? 'n/a'
        : v.toFixed(sigma > 0 ? Math.max(0, Math.min(9, 1 - Math.floor(Math.log10(sigma)))) : fallback);
    const spread = (v) => (v == null ? 'n/a' : v === 0 ? '0' : Number(v.toPrecision(2)).toString());
    const row = (label, key, unit) => [
      `${label} (${unit})`,
      precise(s[key].mean, s[key].sigma),
      spread(s[key].sigma),
    ];
    const day = Math.floor(s.epoch_mjd_utc),
      micros = Math.round((s.epoch_mjd_utc - day) * 86400e6);
    const utc =
      new Date((day - 40587) * 86400e3 + Math.floor(micros / 1000)).toISOString().slice(0, 19) +
      '.' +
      String(micros % 1000000).padStart(6, '0');
    const rows = [
      ['Time (UTC)', utc, `${spread(s.time_sigma_s)}${s.time_sigma_s == null ? '' : ' s'}`],
      row('WGS84 height', 'height_wgs84_km', 'km'),
      row('Latitude', 'latitude_deg', 'deg N'),
      row('Longitude', 'longitude_deg', 'deg E'),
      row('Earth-fixed speed', 'ecef_speed_km_s', 'km/s'),
      row('Flight path angle', 'flight_path_deg', 'deg'),
      row('Heading', 'heading_deg', 'deg'),
      ['1-sigma semi-major axis (m)', format(s.ellipse.semi_major_m), ''],
      ['1-sigma semi-minor axis (m)', format(s.ellipse.semi_minor_m), ''],
      ['Ellipse azimuth (deg)', format(s.ellipse.azimuth_deg), ''],
    ];
    table(
      'entry-statistics',
      'Atmospheric entry',
      ['Parameter', 'Weighted mean / value', '1-sigma spread'],
      rows,
    );
    el('entry-note').textContent =
      `${s.positive_weight_entries} positive-weight entries of ${s.total_entries} archived; effective N = ${format(s.effective_entries)}. Original entry weights, conditional on atmospheric entry; not atmospheric MC resampling. ` +
      'Speed and angles use the Earth-fixed (ECEF) velocity. Flight path angle is relative to the local WGS84 horizontal, negative downward; heading is clockwise from north. ' +
      'Spreads are weighted population standard deviations, not errors of the mean; longitude and heading use wrapped angular residuals. ' +
      'The horizontal ellipse uses local east/north covariance at entry height; azimuth is east of north modulo 180 deg. Its 1-sigma axes are not a 68% joint contour. ' +
      'n/a means undefined or insufficient samples. ' +
      (s.warnings ?? []).join(' ');
  }
  const map = L.map('map', {
    preferCanvas: true,
    zoomControl: true,
    attributionControl: true,
    zoomAnimation: false,
    zoomSnap: 0.1,
    trackResize: false,
  });
  map.setView([0, 0], 2);
  L.control.scale({ imperial: false }).addTo(map);
  map.createPane('cached-imagery').style.zIndex = 200;
  map.createPane('online-imagery').style.zIndex = 201;
  const onlineImagery = data.options.basemap === 'esri_imagery' && !data.options.offline;
  let imageryFailed = false,
    imageryAdded = false;
  const imagery = onlineImagery
    ? L.tileLayer(
        'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}?blankTile=false',
        {
          pane: 'online-imagery',
          maxZoom: 22,
          maxNativeZoom: 19,
          keepBuffer: 2,
          attribution: 'Imagery &copy; Esri, Vantor, Earthstar Geographics, GIS User Community',
        },
      )
        .on('loading', () => {
          imageryFailed = false;
        })
        .on('tileerror', () => {
          imageryFailed = true;
          updateMapMeta();
        })
        .on('load', () => updateMapMeta())
    : null;
  const tileLat = (y) => (Math.atan(Math.sinh(Math.PI * (1 - 2 * y))) * 180) / Math.PI;
  for (const tile of data.tiles) {
    const n = 2 ** tile.z;
    L.imageOverlay(
      tile.url,
      [
        [tileLat((tile.y + 1) / n), (tile.x / n) * 360 - 180],
        [tileLat(tile.y / n), ((tile.x + 1) / n) * 360 - 180],
      ],
      { interactive: false, pane: 'cached-imagery' },
    ).addTo(map);
  }
  if (data.tiles.length) map.attributionControl.addAttribution('Imagery &copy; Esri (embedded tiles)');
  const paletteStops = [
    [68, 1, 84],
    [59, 82, 139],
    [33, 145, 140],
    [94, 201, 98],
    [253, 231, 37],
  ];
  function palette(value) {
    const v = Math.max(0, Math.min(1, value)) * 4,
      i = Math.min(3, Math.floor(v));
    return paletteStops[i].map((c, k) => Math.round(c + (paletteStops[i + 1][k] - c) * (v - i)));
  }
  const paletteTable = Array.from({ length: 256 }, (_, i) => palette(i / 255));
  const grid = data.grid,
    fields = grid ? { count: await R.unpack(grid.count), mass: await R.unpack(grid.mass) } : {};
  if (grid) {
    delete grid.count;
    delete grid.mass;
  }
  const peaks = {};
  for (const q of ['count', 'mass']) {
    peaks[q] = 0;
    for (const v of fields[q] ?? []) peaks[q] = Math.max(peaks[q], v);
  }
  const canvas = L.DomUtil.create('canvas', 'density-canvas');
  map.getPane('overlayPane').append(canvas);
  let scheduled = false;
  function densityImage(ctx, w, h, coordinate, q) {
    const image = ctx.createImageData(w, h),
      field = fields[q],
      peak = peaks[q];
    if (!grid || !field || !(peak > 0)) return image;
    const xs = new Float64Array(w);
    for (let x = 0; x < w; x++) xs[x] = (coordinate(x + 0.5, 0).lng - grid.longitude) / grid.dlon;
    for (let y = 0; y < h; y++) {
      const gy = (coordinate(0, y + 0.5).lat - grid.latitude) / grid.dlat,
        iy = Math.floor(gy),
        fy = gy - iy;
      if (iy < 0 || iy >= grid.rows - 1) continue;
      for (let x = 0; x < w; x++) {
        const gx = xs[x],
          ix = Math.floor(gx),
          fx = gx - ix;
        if (ix < 0 || ix >= grid.columns - 1) continue;
        const i = iy * grid.columns + ix;
        const v =
          ((1 - fy) * ((1 - fx) * field[i] + fx * field[i + 1]) +
            fy * ((1 - fx) * field[i + grid.columns] + fx * field[i + grid.columns + 1])) /
          peak;
        if (v < 2e-4) continue;
        const rgb = paletteTable[Math.min(255, Math.max(0, Math.round(v * 255)))],
          p = (y * w + x) * 4;
        image.data[p] = rgb[0];
        image.data[p + 1] = rgb[1];
        image.data[p + 2] = rgb[2];
        image.data[p + 3] = 255;
      }
    }
    return image;
  }
  function drawDensity() {
    scheduled = false;
    canvas.style.display = grid && el('density').checked ? 'block' : 'none';
    if (!grid || !el('density').checked) return;
    const size = map.getSize(),
      w = Math.round(size.x),
      h = Math.round(size.y);
    if (!w || !h) return;
    canvas.width = w;
    canvas.height = h;
    L.DomUtil.setPosition(canvas, map.containerPointToLayerPoint([0, 0]));
    canvas.style.opacity = el('opacity').value;
    const ctx = canvas.getContext('2d');
    ctx.putImageData(
      densityImage(ctx, w, h, (x, y) => map.containerPointToLatLng([x, y]), el('quantity').value),
      0,
      0,
    );
  }
  function scheduleDensity() {
    if (!scheduled) {
      scheduled = true;
      requestAnimationFrame(drawDensity);
    }
  }
  map.on('moveend zoomend resize', scheduleDensity);
  const coordinatePane = map.createPane('coordinates');
  coordinatePane.style.zIndex = 450;
  coordinatePane.style.pointerEvents = 'none';
  const coordinateCanvas = L.DomUtil.create('canvas', 'coordinate-canvas');
  coordinatePane.append(coordinateCanvas);
  coordinateCanvas.setAttribute('aria-hidden', 'true');
  let gridScheduled = false;
  function drawCoordinates() {
    gridScheduled = false;
    coordinateCanvas.style.display = el('coordinates').checked ? 'block' : 'none';
    if (!el('coordinates').checked || !el('map').clientWidth || !el('map').clientHeight) return;
    const size = map.getSize(),
      w = Math.round(size.x),
      h = Math.round(size.y),
      ratio = Math.min(2, window.devicePixelRatio || 1);
    if (!w || !h) return;
    coordinateCanvas.width = Math.round(w * ratio);
    coordinateCanvas.height = Math.round(h * ratio);
    coordinateCanvas.style.width = w + 'px';
    coordinateCanvas.style.height = h + 'px';
    L.DomUtil.setPosition(coordinateCanvas, map.containerPointToLayerPoint([0, 0]));
    const ctx = coordinateCanvas.getContext('2d');
    ctx.scale(ratio, ratio);
    StrewnGraticule.draw(
      ctx,
      w,
      h,
      (x, y) => map.containerPointToLatLng([x, y]),
      (lat, lon) => {
        const p = map.latLngToContainerPoint([lat, lon]);
        return [p.x, p.y];
      },
    );
  }
  function scheduleCoordinates() {
    if (!gridScheduled) {
      gridScheduled = true;
      requestAnimationFrame(drawCoordinates);
    }
  }
  map.on('move zoomend resize', scheduleCoordinates);
  el('coordinates').onchange = scheduleCoordinates;
  let contours = L.layerGroup().addTo(map),
    points = L.layerGroup().addTo(map);
  const legend = L.control({ position: 'bottomright' });
  legend.onAdd = () => {
    const node = L.DomUtil.create('div', 'map-legend');
    node.id = 'map-legend';
    return node;
  };
  legend.addTo(map);
  function refreshMap() {
    const q = el('quantity').value;
    contours.clearLayers();
    points.clearLayers();
    if (el('hpd').checked)
      L.geoJSON(data.contours, {
        filter: (f) => f.properties.quantity === q,
        style: (f) => ({
          weight: 2,
          color:
            f.properties.hpd > 0.98
              ? '#f75464'
              : f.properties.hpd > 0.7
                ? '#ead53f'
                : data.tiles.length
                  ? '#f4ffff'
                  : '#476670',
          opacity: 1,
        }),
      }).addTo(contours);
    if (el('points').checked)
      for (const p of data.impacts) {
        if (el('point-scope').value === 'selected' && p[5] !== runs[selection]?.run_id) continue;
        const popup = document.createElement('div');
        popup.textContent = `Run ${p[5]}, fragment ${p[6]} | ${format(p[2] * 1000)} g | multiplicity ${format(p[3])} | ${p[0].toFixed(6)}, ${p[1].toFixed(6)}`;
        L.circleMarker([p[0], p[1]], {
          radius: Math.min(7, 2 + Math.log10(1 + p[2] * 1000)),
          weight: 0.7,
          color: '#fff',
          fillColor: '#ca893d',
          fillOpacity: 0.8,
        })
          .bindPopup(popup)
          .addTo(points);
      }
    el('map-legend').innerHTML =
      '<strong>Ensemble probability / cell</strong><div class="gradient"></div><div class="legend-row"><span>0</span><span id="peak"></span></div><div>99% <span style="color:#e44859">&#9473;</span> &nbsp;80% <span style="color:#ba9d00">&#9473;</span> &nbsp;50% <span style="color:#516b75">&#9473;</span> HPD</div>';
    el('peak').textContent = peaks[q].toExponential(3);
    updateMapMeta();
    scheduleDensity();
    scheduleCoordinates();
  }
  function updateMapMeta() {
    const q = el('quantity').value;
    const imageryStatus = onlineImagery
      ? imageryFailed
        ? 'Some online imagery is unavailable; embedded tiles cover the saved field.'
        : 'Online imagery with embedded fallback for the saved field.'
      : data.tiles.length
        ? 'Offline imagery: coverage limited to the embedded tiles.'
        : 'No imagery; coordinates and density are available offline.';
    el('map-meta').textContent = grid
      ? `${q === 'count' ? 'Number' : 'Mass'}-weighted ensemble KDE; cell ${format(grid.cell_m)} m; normalized to 1, conditional on atmospheric entry. ${el('point-scope').value === 'selected' ? 'Point filter: run ' + runs[selection]?.run_id + '. ' : ''}${imageryStatus}`
      : 'No surviving impacts. No probability field is defined.';
    if (grid) {
      const info = grid.metadata[q],
        c = info.kernel_covariance_m2,
        delta = Math.hypot(c[0] - c[2], 2 * c[1]);
      const sigmas = [
        Math.sqrt(Math.max(0, (c[0] + c[2] + delta) / 2)),
        Math.sqrt(Math.max(0, (c[0] + c[2] - delta) / 2)),
      ];
      el('map-meta').textContent +=
        ` ${info.fixed_kernel_fallback ? 'Fixed fallback' : data.options.method} kernel sigma: ${format(sigmas[0])} x ${format(sigmas[1])} m${data.options.method === 'adaptive' && !info.fixed_kernel_fallback ? ' (pilot); local scale ' + data.options.adaptive_limits.join(' - ') : ''}.`;
    }
    const preview = data.impact_preview;
    if (preview && preview.shown < preview.total)
      el('map-meta').textContent +=
        ` Point preview: ${format(preview.shown)} of ${format(preview.total)} tracked impacts, stratified by mass within each run. Density, HPD and distributions use all impacts.`;
  }
  function fit() {
    if (data.bounds && el('map').clientWidth > 36 && el('map').clientHeight > 36)
      map.fitBounds(data.bounds, { padding: [18, 18], animate: false });
  }
  let fitted = false;
  window.addEventListener('resize', () => {
    if (!document.body.classList.contains('print-ready') && active === 'field' && el('map').clientWidth > 36)
      map.invalidateSize({ pan: false });
  });
  el('fit').onclick = fit;
  ['quantity', 'density', 'hpd', 'points', 'point-scope'].forEach((id) => (el(id).onchange = refreshMap));
  el('opacity').oninput = scheduleDensity;
  const mapFlight = StrewnMapFlight.create({
    map,
    L,
    palette,
    choose,
    load: async () => {
      const run = runs[selection],
        f = currentFlight;
      if (!f?.verified) throw new Error(f?.error || 'Trajectory replay unavailable.');
      const values = decoded ?? (await R.unpack(f.samples));
      if (f === currentFlight) decoded = values;
      return {
        run,
        flight: f,
        values,
        stride: data.trace_columns.length,
        threshold: (data.job.physics?.luminous_speed_mps ?? 3000) / 1000,
      };
    },
  });
  async function choose(index) {
    const token = ++selectionRevision;
    mapFlight.pause();
    ++flightRevision;
    selection = Math.max(0, Math.min(runs.length - 1, index));
    decoded = null;
    el('run').value = String(selection);
    el('previous').disabled = selection <= 0;
    el('next').disabled = selection >= runs.length - 1;
    const run = runs[selection];
    el('fragment').replaceChildren(new Option('All fragments', 'all'));
    currentFlight = null;
    const loaded = run ? await R.readFlight(el(run.flight.element).textContent) : null;
    if (token !== selectionRevision) return;
    currentFlight = loaded;
    el('fragment').replaceChildren(
      new Option(currentFlight?.trace_selection?.sampled ? 'All preview tracks' : 'All fragments', 'all'),
    );
    el('map-run').value = String(selection);
    await mapFlight.refresh();
    if (token !== selectionRevision) return;
    for (const t of currentFlight?.traces ?? [])
      el('fragment').add(new Option(`#${t.id} | ${t.end_event} | x${format(t.multiplicity)}`, String(t.id)));
    if (active === 'flight') await flight();
    if (active === 'wind') await wind();
    refreshMap();
  }
  function samples(track, values) {
    const stride = data.trace_columns.length,
      result = [];
    for (let i = 0; i < track.count; i++)
      result.push(Array.from(values.subarray(track.offset + i * stride, track.offset + (i + 1) * stride)));
    result[0] = track.first;
    result[result.length - 1] = track.last;
    return result;
  }
  async function flight() {
    const revision = ++flightRevision;
    const run = runs[selection];
    if (!run) return;
    const f = currentFlight;
    if (!f) return;
    const breakups = f.breakup_rows ?? run.breakups ?? [];
    const description = `Run ${run.run_id} | Entry ${run.entry_id} | Epoch MJD UTC ${run.entry?.epoch_mjd_utc ?? 'n/a'} | Landed ${format(run.mass_balance_kg.landed)} kg | ${format(run.ground_count)} fragments (${run.tracked_count} representatives) | ${run.breakup_count ?? breakups.length} tracked breakups`;
    el('run-meta').textContent = description;
    el('run-print').textContent =
      `Fragments: ${el('fragment').value === 'all' ? (f.trace_selection?.sampled ? 'preview sample' : 'all') : el('fragment').value}; start: ${el('start').value}; time axis: ${el('log-time').checked ? 'log' : 'linear'}.`;
    el('run-error').textContent = f.verified ? '' : `Trajectory replay unavailable: ${f.error}`;
    el('trace-note').textContent = f.verified
      ? `${f.trace_selection?.sampled ? `Displayed ${f.traces.length} of ${f.trace_selection.total} tracks. ` : ''}${data.trace_sampling} Full replay verified against archived impacts, mass budget and breakup events. ${data.range_definition}.`
      : '';
    table(
      'events',
      `Breakup events: run ${run.run_id}`,
      [
        'Fragment',
        'Parent',
        'Generation',
        'Time / s',
        'WGS84 / km',
        'Parent mass / kg',
        'Strength / MPa',
        'Multiplicity',
        'Child groups',
      ],
      [...breakups]
        .sort((a, b) => a.time_s - b.time_s)
        .map((b) => [
          b.fragment_id,
          b.parent_id,
          b.generation,
          format(b.time_s),
          format(b.height_wgs84_m / 1000),
          format(b.mass_kg),
          format(b.strength_pa / 1e6),
          format(b.multiplicity),
          b.groups,
        ]),
    );
    table(
      'parameters',
      'Saved realization parameters',
      ['Parameter', 'Value'],
      Object.entries(run.parameters).map(([k, v]) => [k, format(v)]),
    );
    if (!f.verified) {
      for (const id of ['altitude', 'mass', 'speed', 'range'])
        await plot(id, [], layout('Replay unavailable', 'Time / s', ''));
      return;
    }
    if (!decoded) {
      const values = await R.unpack(f.samples);
      if (revision !== flightRevision) return;
      decoded = values;
    }
    const selected = el('fragment').value,
      tracks = f.traces.filter((t) => selected === 'all' || String(t.id) === selected);
    const first = breakups.length ? breakups.reduce((min, b) => Math.min(min, b.time_s), Infinity) : 0;
    const start = el('start').value === 'breakup' ? first : 0,
      log = el('log-time').checked;
    const paths = tracks.map((t) => ({
      track: t,
      rows: samples(t, decoded).filter((row) => row[0] >= start && (!log || row[0] > 0)),
    }));
    let lo = Infinity,
      hi = -Infinity;
    for (const p of paths)
      for (const row of p.rows)
        if (row[2] > 0) {
          lo = Math.min(lo, Math.log10(row[2]));
          hi = Math.max(hi, Math.log10(row[2]));
        }
    el('mass-scale').textContent = Number.isFinite(lo)
      ? `${format(10 ** lo)} - ${format(10 ** hi)} kg (log)`
      : '';
    const specs = [
      ['altitude', 1, 'Altitude', 'WGS84 height / km'],
      ['mass', 2, 'Fragment mass', 'kg (log)'],
      ['speed', 3, 'Speed', 'ECEF speed / km/s'],
      ['range', 4, 'Ground range', 'Displacement from entry / km'],
    ];
    for (const [id, column, title, unit] of specs) {
      const buckets = Array.from({ length: 16 }, (_, i) => ({
        x: [],
        y: [],
        text: [],
        mode: 'lines',
        type: 'scatter',
        line: { width: 1.2, color: `rgb(${palette(i / 15).join(',')})` },
        hovertemplate: '%{text}<br>%{x:.5g} s<br>%{y:.5g}<extra></extra>',
      }));
      for (const p of paths)
        for (let i = 1; i < p.rows.length; i++) {
          const a = p.rows[i - 1],
            b = p.rows[i],
            mass = Math.sqrt(Math.max(1e-30, a[2] * b[2]));
          const bin = Math.min(
              15,
              Math.max(0, Math.floor(((Math.log10(mass) - lo) / Math.max(1e-9, hi - lo)) * 15)),
            ),
            bucket = buckets[bin];
          bucket.x.push(a[0], b[0], null);
          bucket.y.push(a[column], b[column], null);
          const text = `Fragment ${p.track.id}, x${format(p.track.multiplicity)}`;
          bucket.text.push(text, text, null);
        }
      const events = f.breakups.filter(
        (b) =>
          (selected === 'all' || String(b.fragment_id) === selected) &&
          b.state[0] >= start &&
          (!log || b.state[0] > 0),
      );
      buckets.push({
        x: events.map((b) => b.state[0]),
        y: events.map((b) => b.state[column]),
        mode: 'markers',
        marker: { symbol: 'x', size: 7, color: '#182c35' },
        text: events.map((b) => `Breakup #${b.fragment_id}`),
        hovertemplate: '%{text}<br>%{x:.6g} s<br>%{y:.6g}<extra></extra>',
      });
      const options = layout(title, `Time since entry / s${log ? ' (log)' : ''}`, unit);
      options.xaxis.type = log ? 'log' : 'linear';
      if (column === 2) options.yaxis.type = 'log';
      if (revision !== flightRevision) return;
      await plot(
        id,
        buckets.filter((b) => b.x.length),
        options,
      );
    }
  }
  async function wind() {
    const run = runs[selection];
    if (!run || !currentFlight) return;
    const rows = currentFlight.wind ?? [];
    el('wind-meta').textContent =
      `Run ${run.run_id}, entry ${run.entry_id}. ENU wind at the entry subpoint; multiplier ${format(run.parameters.wind_scale)} applied. Direction TO, not meteorological FROM.${run.flight.verified ? '' : ' Replay unavailable.'}`;
    const traces = ['East (U)', 'North (V)', 'Up (W)', 'Horizontal speed'].map((name, i) => ({
      x: rows.map((r) => (i === 3 ? Math.hypot(r[1], r[2]) : r[i + 1])),
      y: rows.map((r) => r[0]),
      name,
      mode: 'lines',
      line: { color: colors[i] },
    }));
    await plot(
      'wind-components',
      traces,
      layout('Wind components', 'm/s', 'WGS84 height / km', {
        showlegend: true,
        legend: { orientation: 'h', y: 1.02, yanchor: 'bottom' },
        margin: { l: 75, r: 25, t: 100, b: 65 },
      }),
    );
    const profile = R.windProfile(rows),
      top = Math.max(1, ...profile.map((p) => p.height)),
      maxSpeed = Math.max(1, ...profile.map((p) => p.speed));
    const buckets = Array.from({ length: 32 }, (_, i) => ({
      type: 'scatterpolar',
      mode: 'lines',
      r: [],
      theta: [],
      showlegend: false,
      hoverinfo: 'skip',
      line: { width: 2.5, color: `rgb(${palette(i / 31).join(',')})` },
    }));
    for (let i = 1; i < profile.length; i++) {
      const a = profile[i - 1],
        b = profile[i];
      if (a.direction === null || b.direction === null) continue;
      const trace = buckets[Math.min(31, Math.floor(((a.speed + b.speed) / 2 / maxSpeed) * 31))];
      trace.r.push(a.height, b.height, null);
      trace.theta.push(a.direction, b.direction, null);
    }
    const polarTraces = buckets.filter((b) => b.r.length);
    polarTraces.push({
      type: 'scatterpolar',
      mode: 'markers',
      r: profile.map((p) => p.height),
      theta: profile.map((p) => p.direction),
      showlegend: false,
      customdata: profile.map((p) => p.speed),
      marker: {
        size: 4,
        color: profile.map((p) => p.speed),
        cmin: 0,
        cmax: maxSpeed,
        colorscale: paletteStops.map((rgb, i) => [i / 4, `rgb(${rgb.join(',')})`]),
        showscale: true,
        colorbar: { title: { text: 'Wind speed (m/s)', side: 'right' }, thickness: 14, len: 0.85, x: 1.02 },
      },
      hovertemplate:
        'Height %{r:.2f} km<br>Direction TO %{theta:.1f} deg<br>Speed %{customdata:.2f} m/s<extra></extra>',
    });
    for (const [point, name, symbol] of [
      [profile[0], 'Surface', 'square'],
      [profile.at(-1), 'Top', 'triangle-up'],
    ])
      if (point) {
        const label = name === 'Surface' && point.height !== 0 ? 'Bottom' : name;
        polarTraces.push({
          type: 'scatterpolar',
          mode: 'markers',
          r: [point.height],
          theta: [point.direction ?? (point.height === 0 ? 0 : null)],
          name: `${label} (${format(point.height)} km)`,
          showlegend: true,
          marker: { symbol, size: 10, color: '#182c35' },
          hovertemplate: `${label}<br>Height %{r:.2f} km<extra></extra>`,
        });
      }
    await plot('wind-direction', polarTraces, {
      title: { text: 'Wind direction vs altitude (km)', x: 0.08, font: { size: 15 } },
      font: { family: 'system-ui, Segoe UI, sans-serif', size: 12, color: '#334a54' },
      paper_bgcolor: '#fff',
      margin: { l: 55, r: 100, t: 70, b: 75 },
      showlegend: true,
      legend: { orientation: 'h', x: 0, y: -0.13 },
      polar: {
        bgcolor: '#fff',
        radialaxis: { range: [0, top * 1.08], dtick: top > 20 ? 10 : 5, angle: 45, gridcolor: '#d7e0e4' },
        angularaxis: {
          rotation: 90,
          direction: 'clockwise',
          tickvals: [0, 45, 90, 135, 180, 225, 270, 315],
          ticktext: ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'],
          gridcolor: '#d7e0e4',
        },
      },
    });
  }
  async function view(name) {
    if (name !== 'field') mapFlight.pause();
    active = name;
    document.querySelectorAll('.view').forEach((n) => n.classList.toggle('active', n.id === name));
    document
      .querySelectorAll('nav button')
      .forEach((n) => n.setAttribute('aria-selected', String(n.dataset.view === name)));
    if (name === 'field') {
      map.invalidateSize({ pan: false });
      if (!fitted) {
        fit();
        fitted = true;
      }
      if (imagery && !imageryAdded) {
        imagery.addTo(map);
        imageryAdded = true;
      }
      refreshMap();
      mapFlight.draw();
    }
    if (name === 'flight') await flight();
    if (name === 'wind') await wind();
    for (const id of plots) if (el(id).offsetWidth) await Plotly.Plots.resize(el(id));
  }
  document.querySelectorAll('nav button').forEach((b) => (b.onclick = () => view(b.dataset.view)));
  el('run').onchange = () => choose(Number(el('run').value));
  el('previous').onclick = () => choose(selection - 1);
  el('next').onclick = () => choose(selection + 1);
  ['fragment', 'start', 'log-time'].forEach((id) => (el(id).onchange = flight));
  function printError(error) {
    el('print-error').textContent = 'Print failed: ' + error.message;
  }
  function finishPrint() {
    document.body.classList.remove('print-ready');
    el('print-pages').replaceChildren();
    el('print').disabled = false;
  }
  window.addEventListener('afterprint', finishPrint);
  el('print').onclick = async () => {
    if (el('print').disabled) return;
    mapFlight.pause();
    el('print').disabled = true;
    el('print-error').textContent = '';
    try {
      const scope = el('print-scope').value === 'all' ? ['overview', 'field', 'flight', 'wind'] : [active];
      if (scope.includes('flight')) await flight();
      if (scope.includes('wind')) await wind();
      await StrewnPrint.prepare({
        data,
        scope,
        plots: printPlots,
        densityImage,
        palette,
        peaks,
        runId: runs[selection]?.run_id,
      });
      document.body.classList.add('print-ready');
      await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
      window.print();
    } catch (error) {
      printError(error);
      finishPrint();
    } finally {
      el('print').disabled = false;
    }
  };
  window.addEventListener('keydown', (event) => {
    if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === 'p') {
      event.preventDefault();
      el('print').onclick();
    }
  });
  el('subtitle').textContent =
    `${data.summary.completed} / ${data.summary.requested} realizations | ${(data.impact_preview?.total ?? data.impacts.length).toLocaleString('en-US')} tracked impacts | ${data.summary.converged ? 'JS converged' : 'Run limit reached'} | ${format(data.summary.elapsed_s)} s`;
  el('warnings').textContent = data.warnings.join('\n');
  el('configuration').textContent = JSON.stringify(
    { impacts_fnv1a64: data.impacts_fnv1a64, job: data.job, map_options: data.options },
    null,
    2,
  );
  el('footer').textContent =
    `Numerical archive: impacts FNV-1a ${data.impacts_fnv1a64}. Trace preview: float32; exact event endpoints retained. Charts: Plotly.js 3.1.0 (MIT); maps: Leaflet (BSD-2-Clause).`;
  runs.forEach((r, i) =>
    el('run').add(
      new Option(
        `Run ${r.run_id} | ${format(r.mass_balance_kg.landed)} kg | N=${format(r.ground_count)}`,
        String(i),
      ),
    ),
  );
  runs.forEach((r, i) =>
    el('map-run').add(new Option(`Run ${r.run_id} | ${format(r.mass_balance_kg.landed)} kg`, String(i))),
  );
  if (!runs.length) {
    el('run-error').textContent = 'Realization archive unavailable.';
    el('run').disabled = true;
    el('previous').disabled = true;
    el('next').disabled = true;
  }
  await overview();
  await choose(0);
  await flight();
  await wind();
  refreshMap();
  document.body.dataset.ready = 'true';
})().catch((error) => {
  document.getElementById('warnings').textContent += '\nReport error: ' + error.message;
  console.error(error);
});
