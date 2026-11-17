(() => {
  'use strict';
  const data = JSON.parse(document.getElementById('report-data').textContent);
  const names = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Moon'];
  const el = (id) => document.getElementById(id);
  const svg = el('plot'),
    tip = el('tooltip'),
    body = el('body'),
    units = el('units');
  const NS = 'http://www.w3.org/2000/svg',
    samples = data.distanceSamples || [];
  const pad = (n) => String(n).padStart(2, '0');
  const date = (mjd) => {
    const d = new Date((mjd - 40587) * 86400000);
    return `${d.getUTCFullYear()}-${pad(d.getUTCMonth() + 1)}-${pad(d.getUTCDate())} ${pad(d.getUTCHours())}:${pad(d.getUTCMinutes())}:${pad(d.getUTCSeconds())}`;
  };
  const factor = () => (units.value === 'km' ? 149597870.7 : units.value === 'ld' ? 149597870.7 / 384400 : 1);
  const unitLabel = () => (units.value === 'ld' ? 'LD' : units.value);
  const num = (n) => {
    if (n == null || !Number.isFinite(n))
      return '—';
    const value = n * factor();
    return units.value === 'km' && value >= 1000
      ? value.toLocaleString('en-US', { maximumFractionDigits: 3 })
      : value.toPrecision(10);
  };
  const speed = (n) => (n == null || !Number.isFinite(n) ? '—' : n.toFixed(5));
  const sigma = (d) =>
    d == null || !Number.isFinite(d)
      ? ''
      : ` ± ${d * 86400 < 60 ? (d * 86400).toFixed(1) + ' s' : (d * 1440).toFixed(1) + ' min'}`;
  const node = (tag, attrs, parent = svg) => {
    const e = document.createElementNS(NS, tag);
    Object.entries(attrs).forEach(([k, v]) => e.setAttribute(k, v));
    parent.appendChild(e);
    return e;
  };
  const text = (x, y, value, anchor = 'middle') => {
    const t = node('text', { x, y, 'text-anchor': anchor });
    t.textContent = value;
    return t;
  };
  names.forEach((name, index) => {
    const option = document.createElement('option');
    option.value = index;
    option.textContent = name;
    body.appendChild(option);
  });
  body.value = '2';
  el('object').textContent = data.object || 'Close Approaches';
  el('interval').textContent = `${date(data.startMjdTdb)} — ${date(data.endMjdTdb)} TDB`;
  let log = true,
    selected = -1,
    left = data.startMjdTdb,
    right = data.endMjdTdb;
  const plot = { left: 80, right: 1180, top: 15, bottom: 295 };
  function resizeValueLabels() {
    const width = svg.getBoundingClientRect().width;
    if (width > 0)
      svg.style.setProperty('--value-label-size', `${10 * Math.min(1, 1200 / width)}px`);
  }
  function drawTable() {
    const tbody = el('rows');
    tbody.replaceChildren();
    let count = 0;
    data.events.forEach((event, index) => {
      if (event.body !== names[+body.value])
        return;
      count++;
      const tr = document.createElement('tr');
      tr.tabIndex = 0;
      if (index === selected)
        tr.className = 'selected';
      [
        date(event.mjdTdb) + sigma(event.time3SigmaDays),
        event.body,
        num(event.distanceAu),
        num(event.minimum3SigmaAu),
        num(event.maximum3SigmaAu),
        speed(event.relativeSpeedKmS),
        speed(event.vInfinityKmS),
      ].forEach((value) => {
        const td = document.createElement('td');
        td.textContent = value;
        tr.appendChild(td);
      });
      tr.title = event.uncertaintyNote || '';
      tr.addEventListener('click', () => {
        selected = index;
        draw();
      });
      tr.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') {
          e.preventDefault();
          selected = index;
          draw();
        }
      });
      tbody.appendChild(tr);
    });
    el('count').textContent = `${count} events`;
  }
  function draw() {
    const idx = +body.value;
    el('chart-title').textContent = `Distance to ${names[idx]} (${unitLabel()})`;
    el('log').classList.toggle('active', log);
    el('linear').classList.toggle('active', !log);
    el('hill-legend').hidden = !el('hill').checked;
    svg.replaceChildren();
    drawTable();
    if (!samples.length) {
      svg.setAttribute('hidden', '');
      el('empty').hidden = false;
      el('empty').textContent = data.stoppedAtSurface
        ? 'Distance curve unavailable for a surface crossing.'
        : 'No distance samples available.';
      return;
    }
    svg.removeAttribute('hidden');
    el('empty').hidden = true;
    let first = 0,
      last = samples.length - 1;
    while (first < samples.length - 1 && samples[first + 1][0] < left)
      first++;
    while (last > 0 && samples[last - 1][0] > right)
      last--;
    const visible = samples.slice(first, last + 1);
    let values = visible.map((s) => s[1 + idx]).filter((x) => Number.isFinite(x) && x > 0);
    if (el('hill').checked)
      values = values.concat(visible.map((s) => s[10 + idx]).filter((x) => Number.isFinite(x) && x > 0));
    if (!values.length)
      return;
    let min = Infinity,
      max = -Infinity;
    values.forEach((value) => {
      min = Math.min(min, value);
      max = Math.max(max, value);
    });
    if (log) {
      min = Math.log10(min);
      max = Math.log10(max);
      const margin = Math.max((max - min) * 0.08, 0.08);
      min -= margin;
      max += margin;
    } else {
      min = 0;
      max = Math.max(max * 1.08, 1e-12);
    }
    const X = (t) => plot.left + ((t - left) / (right - left)) * (plot.right - plot.left);
    const Y = (v) =>
      plot.bottom -
      (((log ? Math.log10(Math.max(v, 1e-30)) : v) - min) / (max - min)) * (plot.bottom - plot.top);
    for (let i = 0; i <= 5; i++) {
      const x = plot.left + ((plot.right - plot.left) * i) / 5,
        t = left + ((right - left) * i) / 5;
      node('line', { x1: x, y1: plot.top, x2: x, y2: plot.bottom, class: 'grid' });
      text(x, 317, date(t).slice(0, 10), i === 5 ? 'end' : 'middle').setAttribute('class', 'date-label');
    }
    for (let i = 0; i <= 5; i++) {
      const y = plot.bottom - ((plot.bottom - plot.top) * i) / 5,
        raw = min + ((max - min) * i) / 5,
        v = (log ? 10 ** raw : raw) * factor();
      node('line', { x1: plot.left, y1: y, x2: plot.right, y2: y, class: 'grid' });
      text(70, y, v.toPrecision(3), 'end').setAttribute('class', 'value-label');
    }
    node('line', { x1: plot.left, y1: plot.bottom, x2: plot.right, y2: plot.bottom, class: 'axis' });
    const pathFor = (col) =>
      visible.map((s, i) => `${i ? 'L' : 'M'}${X(s[0]).toFixed(2)} ${Y(s[col]).toFixed(2)}`).join(' ');
    const defs = node('defs', {}),
      clip = node('clipPath', { id: 'clip' }, defs);
    node(
      'rect',
      { x: plot.left, y: plot.top, width: plot.right - plot.left, height: plot.bottom - plot.top },
      clip,
    );
    const group = node('g', { 'clip-path': 'url(#clip)' });
    node('path', { d: pathFor(1 + idx), class: 'trajectory' }, group);
    if (el('hill').checked)
      node('path', { d: pathFor(10 + idx), class: 'hill' }, group);
    data.events.forEach((event, index) => {
      if (event.body !== names[idx] || event.mjdTdb < left || event.mjdTdb > right || !(event.distanceAu > 0))
        return;
      const c = node(
        'circle',
        {
          cx: X(event.mjdTdb),
          cy: Y(event.distanceAu),
          r: index === selected ? 5 : 3.5,
          class: index === selected ? 'event selected' : 'event',
        },
        group,
      );
      const title = node('title', {}, c);
      title.textContent = `${date(event.mjdTdb)} TDB · ${num(event.distanceAu)} ${unitLabel()}`;
      c.addEventListener('click', () => {
        selected = index;
        draw();
        const row = document.querySelector('tr.selected');
        if (row)
          row.scrollIntoView({ block: 'nearest' });
      });
    });
  }
  function zoom(center, multiplier) {
    const span = Math.min(
      data.endMjdTdb - data.startMjdTdb,
      Math.max(1 / 86400, (right - left) * multiplier),
    );
    const ratio = (center - left) / (right - left);
    left = center - span * ratio;
    right = left + span;
    if (left < data.startMjdTdb) {
      left = data.startMjdTdb;
      right = left + span;
    }
    if (right > data.endMjdTdb) {
      right = data.endMjdTdb;
      left = right - span;
    }
    draw();
  }
  svg.addEventListener(
    'wheel',
    (e) => {
      e.preventDefault();
      const r = svg.getBoundingClientRect(),
        x = ((e.clientX - r.left) / r.width) * 1200,
        t = left + Math.max(0, Math.min(1, (x - plot.left) / (plot.right - plot.left))) * (right - left);
      zoom(t, e.deltaY < 0 ? 0.7 : 1.4);
    },
    { passive: false },
  );
  svg.addEventListener('mousemove', (e) => {
    if (!samples.length)
      return;
    const r = svg.getBoundingClientRect(),
      x = ((e.clientX - r.left) / r.width) * 1200;
    if (x < plot.left || x > plot.right) {
      tip.style.display = 'none';
      return;
    }
    const t = left + ((x - plot.left) / (plot.right - plot.left)) * (right - left);
    let lo = 0,
      hi = samples.length - 1;
    while (lo < hi) {
      const m = (lo + hi) >> 1;
      if (samples[m][0] < t)
        lo = m + 1;
      else
        hi = m;
    }
    if (lo > 0 && Math.abs(samples[lo - 1][0] - t) < Math.abs(samples[lo][0] - t))
      lo--;
    tip.textContent = `${date(samples[lo][0])} TDB  ·  ${num(samples[lo][1 + +body.value])} ${unitLabel()}`;
    tip.style.left = Math.min(e.clientX - r.left + 12, r.width - 210) + 'px';
    tip.style.top = '10px';
    tip.style.display = 'block';
  });
  svg.addEventListener('mouseleave', () => (tip.style.display = 'none'));
  body.addEventListener('change', () => {
    selected = -1;
    draw();
  });
  units.addEventListener('change', draw);
  el('hill').addEventListener('change', draw);
  el('log').addEventListener('click', () => {
    log = true;
    draw();
  });
  el('linear').addEventListener('click', () => {
    log = false;
    draw();
  });
  el('reset').addEventListener('click', () => {
    left = data.startMjdTdb;
    right = data.endMjdTdb;
    draw();
  });
  window.addEventListener('resize', resizeValueLabels);
  resizeValueLabels();
  draw();
})();
