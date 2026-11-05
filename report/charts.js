// Gráficas de la página, dibujadas en SVG a partir de los datos de results/fig_*.json (report/build.py los incrusta).
// Colores por variables CSS (--c1, --c2, --c3, --ink, --muted, --rule, --sheet), así que siguen el modo claro u oscuro.
(() => {
  const DATA = JSON.parse(document.getElementById('figdata').textContent);
  const NS = 'http://www.w3.org/2000/svg';
  const PI = Math.PI;

  const svgEl = (tag, attrs = {}, parent) => {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  };
  const htmlEl = (tag, cls, parent) => {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (parent) parent.appendChild(e);
    return e;
  };

  // "m_{eff}^{2}" -> texto con subíndices y superíndices, en SVG o en HTML
  const PARTS = /([_^])\{([^}]*)\}/g;
  function richSvg(t, s) {
    let i = 0, m;
    PARTS.lastIndex = 0;
    while ((m = PARTS.exec(s))) {
      if (m.index > i) t.appendChild(document.createTextNode(s.slice(i, m.index)));
      const ts = svgEl('tspan', { 'baseline-shift': m[1] === '_' ? 'sub' : 'super', 'font-size': '75%' }, t);
      ts.textContent = m[2];
      i = PARTS.lastIndex;
    }
    if (i < s.length) t.appendChild(document.createTextNode(s.slice(i)));
    return t;
  }
  function richHtml(e, s) {
    let i = 0, m;
    PARTS.lastIndex = 0;
    while ((m = PARTS.exec(s))) {
      if (m.index > i) e.appendChild(document.createTextNode(s.slice(i, m.index)));
      htmlEl(m[1] === '_' ? 'sub' : 'sup', '', e).textContent = m[2];
      i = PARTS.lastIndex;
    }
    if (i < s.length) e.appendChild(document.createTextNode(s.slice(i)));
    return e;
  }

  const minus = s => s.replace(/-/g, '−');
  const fmt = (v, d = 3) => minus(Number(v).toFixed(d));
  const sci = v => {
    const e = Math.floor(Math.log10(v));
    const sup = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
    return e === 0 ? '1' : '10' + String(e).split('').map(c => sup[c]).join('');
  };

  function niceTicks(a, b, n) {
    const step0 = (b - a) / n, p = 10 ** Math.floor(Math.log10(step0)), f = step0 / p;
    const step = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
    const out = [];
    for (let v = Math.ceil(a / step) * step; v <= b + step * 1e-9; v += step) out.push(+v.toFixed(10));
    return out.map(v => [v, minus(String(+v.toPrecision(6)))]);
  }

  function scale(ax, r0, r1) {
    if (ax.log) {
      const l0 = Math.log10(ax.min), l1 = Math.log10(ax.max);
      const f = v => r0 + (Math.log10(Math.max(v, ax.min * 1e-3)) - l0) / (l1 - l0) * (r1 - r0);
      f.inv = p => 10 ** (l0 + (p - r0) / (r1 - r0) * (l1 - l0));
      return f;
    }
    const f = v => r0 + (v - ax.min) / (ax.max - ax.min) * (r1 - r0);
    f.inv = p => ax.min + (p - r0) / (r1 - r0) * (ax.max - ax.min);
    return f;
  }

  const cssVar = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  function rgb(hex) {
    const c = document.createElement('canvas').getContext('2d');
    c.fillStyle = hex;
    const h = c.fillStyle;   // normaliza a #rrggbb
    return [1, 3, 5].map(i => parseInt(h.slice(i, i + 2), 16));
  }

  // clave de leyenda: línea, línea discontinua, área, punto o anillo
  function key(kind, color, parent) {
    const s = svgEl('svg', { viewBox: '0 0 18 10', width: 18, height: 10, 'aria-hidden': 'true' }, parent);
    const st = `stroke:var(${color});fill:var(${color})`;
    if (kind === 'area') svgEl('rect', { x: 1, y: 1, width: 16, height: 8, rx: 2, style: `fill:var(${color});opacity:.25` }, s);
    else if (kind === 'dot') svgEl('circle', { cx: 9, cy: 5, r: 4, style: st }, s);
    else if (kind === 'ring') svgEl('circle', { cx: 9, cy: 5, r: 3.5, style: `stroke:var(${color});fill:var(--sheet);stroke-width:2` }, s);
    else if (kind === 'heat') svgEl('rect', { x: 1, y: 1, width: 16, height: 8, rx: 2, style: `fill:url(#heatkey)` }, s);
    else svgEl('line', { x1: 1, y1: 5, x2: 17, y2: 5, style: `stroke:var(${color});stroke-width:2;stroke-linecap:round` + (kind === 'dash' ? ';stroke-dasharray:4 3' : '') }, s);
    return s;
  }

  // ------------------------------------------------------------------ motor
  function render(box, spec) {
    box.textContent = '';
    const head = htmlEl('div', 'chead', box);
    richHtml(htmlEl('h4', '', head), spec.title);
    if (spec.sub) richHtml(htmlEl('p', 'sub', head), spec.sub);
    if (spec.legend) {
      const lg = htmlEl('div', 'legend', box);
      for (const [kind, color, name] of spec.legend) {
        const it = htmlEl('span', '', lg);
        if (kind === 'heat') {
          const bar = htmlEl('i', 'heatbar', it);
          bar.style.background = `linear-gradient(90deg, var(--c1), var(--sheet))`;
        } else key(kind, color, it);
        richHtml(htmlEl('span', '', it), name);
      }
    }

    const W = Math.max(260, Math.round(box.clientWidth - 2)), H = spec.h || 290;
    const m = { l: spec.ml || 52, r: 14, t: 12, b: spec.mb || 42 };
    const svg = svgEl('svg', { class: 'plot', viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': spec.aria }, box);
    const X = scale(spec.x, m.l, W - m.r), Y = scale(spec.y, H - m.b, m.t);
    const id = 'c' + Math.random().toString(36).slice(2, 8);
    const defs = svgEl('defs', {}, svg);
    const cp = svgEl('clipPath', { id: id + 'clip' }, defs);
    svgEl('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b }, cp);

    // rejilla y ejes
    const g = svgEl('g', { class: 'grid' }, svg);
    const xt = spec.x.ticks || niceTicks(spec.x.min, spec.x.max, Math.max(3, Math.floor((W - m.l) / 80)));
    const yt = spec.y.ticks || niceTicks(spec.y.min, spec.y.max, 5);
    for (const [v, lab] of yt) {
      const y = Y(v);
      svgEl('line', { x1: m.l, x2: W - m.r, y1: y, y2: y }, g);
      svgEl('text', { x: m.l - 7, y: y + 3.5, 'text-anchor': 'end', class: 'tick' }, svg).textContent = lab;
    }
    for (const [v, lab] of xt) {
      const x = X(v);
      if (!spec.x.nogrid) svgEl('line', { x1: x, x2: x, y1: m.t, y2: H - m.b }, g);
      richSvg(svgEl('text', { x, y: H - m.b + 15, 'text-anchor': 'middle', class: 'tick' }, svg), lab);
    }
    svgEl('line', { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: 'axis' }, svg);
    richSvg(svgEl('text', { x: (m.l + W - m.r) / 2, y: H - 6, 'text-anchor': 'middle', class: 'alabel' }, svg), spec.x.label);
    richSvg(svgEl('text', { x: 12, y: (m.t + H - m.b) / 2, 'text-anchor': 'middle', class: 'alabel',
      transform: `rotate(-90 12 ${(m.t + H - m.b) / 2})` }, svg), spec.y.label);

    const plot = svgEl('g', { 'clip-path': `url(#${id}clip)` }, svg);
    const over = svgEl('g', {}, svg);
    const pts = [], cross = [];
    const path = d => d.map(([x, y], i) => (i ? 'L' : 'M') + X(x).toFixed(1) + ' ' + Y(y).toFixed(1)).join('');

    for (const L of spec.layers) {
      const col = `var(${L.color || '--ink'})`;
      if (L.type === 'heat') {
        const cv = document.createElement('canvas');
        cv.width = L.nx; cv.height = L.ny;
        const ctx = cv.getContext('2d'), img = ctx.createImageData(L.nx, L.ny);
        const a = rgb(cssVar('--c1')), b = rgb(cssVar('--sheet'));
        for (let j = 0; j < L.ny; j++) for (let i = 0; i < L.nx; i++) {
          const t = Math.min(1, Math.max(0, (L.f[j * L.nx + i] / 10 - L.vmin) / (L.vmax - L.vmin))) ** 1.7;
          const o = ((L.ny - 1 - j) * L.nx + i) * 4;
          for (let c = 0; c < 3; c++) img.data[o + c] = a[c] + (b[c] - a[c]) * t;
          img.data[o + 3] = 255;
        }
        ctx.putImageData(img, 0, 0);
        svgEl('image', { href: cv.toDataURL(), x: X(L.x[0]), y: Y(L.y[1]), width: X(L.x[1]) - X(L.x[0]),
          height: Y(L.y[0]) - Y(L.y[1]), preserveAspectRatio: 'none', style: 'image-rendering:auto' }, plot);
        spec.heat = L;
      } else if (L.type === 'area') {
        const up = L.data.map(p => [p[0], p[2]]), lo = L.data.map(p => [p[0], p[1]]).reverse();
        svgEl('path', { d: path(up) + path(lo).replace('M', 'L') + 'Z', style: `fill:${col};opacity:${L.opacity || 0.14}` }, plot);
      } else if (L.type === 'line') {
        svgEl('path', { d: path(L.data), class: 'line', style: `stroke:${col};stroke-width:${L.width || 2}` + (L.dash ? `;stroke-dasharray:${L.dash}` : '') }, plot);
        if (L.name) cross.push(L);
      } else if (L.type === 'hline' || L.type === 'vline') {
        const h = L.type === 'hline';
        const a = h ? { x1: m.l, x2: W - m.r, y1: Y(L.at), y2: Y(L.at) } : { x1: X(L.at), x2: X(L.at), y1: m.t, y2: H - m.b };
        svgEl('line', { ...a, style: `stroke:${col};stroke-width:${L.width || 1.2}` + (L.dash ? `;stroke-dasharray:${L.dash}` : '') }, plot);
        if (L.label) {
          const t = h ? { x: L.lx != null ? X(L.lx) : W - m.r - 4, y: Y(L.at) - 5, 'text-anchor': L.lx != null ? 'start' : 'end' }
                      : { x: X(L.at) + (L.left ? -5 : 5), y: L.ly != null ? Y(L.ly) : m.t + 11, 'text-anchor': L.left ? 'end' : 'start' };
          richSvg(svgEl('text', { ...t, class: 'dlabel' }, over), L.label);
        }
      } else if (L.type === 'points') {
        for (const p of L.data) {
          const cx = X(p[0]), cy = Y(p[1]);
          if (L.err) svgEl('line', { x1: cx, x2: cx, y1: Y(p[1] - p[2]), y2: Y(p[1] + p[2]), style: `stroke:${col};stroke-width:1.5;opacity:.55` }, plot);
          svgEl('circle', { cx, cy, r: L.r || 4, class: L.hollow ? 'ring' : 'dot', style: L.hollow ? `stroke:${col}` : `fill:${col}` }, over);
          if (L.tip) pts.push({ cx, cy, color: L.color, tip: L.tip(p) });
        }
      } else if (L.type === 'text') {
        richSvg(svgEl('text', { x: X(L.x) + (L.dx || 0), y: Y(L.y) + (L.dy || 0), 'text-anchor': L.anchor || 'start', class: 'dlabel' }, over), L.label);
      }
    }

    // capa de interacción: punto más cercano (radio 18 px), si no, cruz vertical o valor del mapa
    const tip = htmlEl('div', 'tip', box);
    tip.hidden = true;
    const hair = svgEl('line', { y1: m.t, y2: H - m.b, class: 'hair', visibility: 'hidden' }, svg);
    const hit = svgEl('rect', { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: 'transparent' }, svg);
    const ring = svgEl('circle', { r: 7, class: 'hover', visibility: 'hidden' }, svg);

    function show(rows, px, py) {
      tip.textContent = '';
      for (const [color, value, name] of rows) {
        const row = htmlEl('div', 'row', tip);
        if (color) htmlEl('i', 'k', row).style.background = `var(${color})`;
        richHtml(htmlEl('b', '', row), value);
        if (name) richHtml(htmlEl('span', 'nm', row), name);
      }
      tip.hidden = false;
      const bw = box.clientWidth, tw = tip.offsetWidth, th = tip.offsetHeight;
      const svgTop = svg.getBoundingClientRect().top - box.getBoundingClientRect().top;
      const k = svg.getBoundingClientRect().width / W;
      let left = px * k + 14, top = svgTop + py * k - th - 10;
      if (left + tw > bw - 4) left = px * k - tw - 14;
      if (top < 4) top = svgTop + py * k + 14;
      tip.style.left = Math.max(4, left) + 'px';
      tip.style.top = top + 'px';
    }
    function hide() { tip.hidden = true; hair.setAttribute('visibility', 'hidden'); ring.setAttribute('visibility', 'hidden'); }

    function move(ev) {
      const r = svg.getBoundingClientRect(), k = W / r.width;
      const px = (ev.clientX - r.left) * k, py = (ev.clientY - r.top) * k;
      let best = null, bd = 18 * k;
      for (const p of pts) { const d = Math.hypot(p.cx - px, p.cy - py); if (d < bd) { bd = d; best = p; } }
      if (best) {
        hair.setAttribute('visibility', 'hidden');
        Object.entries({ cx: best.cx, cy: best.cy, visibility: 'visible' }).forEach(([a, v]) => ring.setAttribute(a, v));
        ring.style.stroke = `var(${best.color || '--ink'})`;
        return show(best.tip, best.cx, best.cy);
      }
      ring.setAttribute('visibility', 'hidden');
      const xv = X.inv(px);
      if (cross.length) {
        const rows = [];
        let snap = null;
        for (const L of cross) {
          let bi = 0;
          for (let i = 1; i < L.data.length; i++) if (Math.abs(L.data[i][0] - xv) < Math.abs(L.data[bi][0] - xv)) bi = i;
          snap = snap ?? L.data[bi][0];
          rows.push([L.color, L.fmt ? L.fmt(L.data[bi][1]) : fmt(L.data[bi][1]), L.name]);
        }
        hair.setAttribute('x1', X(snap)); hair.setAttribute('x2', X(snap)); hair.setAttribute('visibility', 'visible');
        return show([[null, spec.x.fmt ? spec.x.fmt(snap) : fmt(snap), spec.x.name || '']].concat(rows), X(snap), py);
      }
      if (spec.heat) {
        const L = spec.heat, yv = Y.inv(py);
        const i = Math.round((xv - L.x[0]) / (L.x[1] - L.x[0]) * (L.nx - 1)), j = Math.round((yv - L.y[0]) / (L.y[1] - L.y[0]) * (L.ny - 1));
        if (i < 0 || j < 0 || i >= L.nx || j >= L.ny) return hide();
        return show([[null, `ω̄ = ${fmt(xv, 2)} ${yv < 0 ? '−' : '+'} ${Math.abs(yv).toFixed(2)}i`, ''],
          [null, fmt(L.f[j * L.nx + i] / 10, 1), 'log_{10}|f|']], px, py);
      }
      hide();
    }
    hit.addEventListener('pointermove', move);
    hit.addEventListener('pointerdown', move);
    hit.addEventListener('pointerleave', hide);
  }

  // ------------------------------------------------------------------ las figuras
  const zip = (a, b) => a.map((v, i) => [v, b[i]]);
  const SPECS = {};

  // masa efectiva frente a theta
  {
    const d = DATA.meff, v10 = d.a10[0], vpi = d.a10[d.a10.length - 1];
    SPECS.meff = {
      title: 'Masa efectiva del axión con las funciones completas', sub: 'β = 0.88; el modelo simplificado usa el valor en θ = 0',
      aria: `m_eff² / m_a² frente a θ. Para α = 10 vale ${fmt(v10, 0)} en θ = 0 y +${vpi.toFixed(0)} en θ = π.`,
      legend: [['line', '--c1', 'α = n_{b}/n_{c} = 10 (interior)'], ['line', '--c2', 'α = 0.5 (fuera de la capa)']],
      x: { min: -0.12, max: PI + 0.12, ticks: [[0, '0'], [PI / 2, 'π/2'], [PI, 'π']], label: 'θ = a/f_{a}  (fondo en el que se linealiza)', fmt: v => 'θ = ' + (v / PI).toFixed(2) + 'π' },
      y: { min: -12, max: 30, label: 'm_{eff}^{2} / m_{a}^{2}' },
      layers: [
        { type: 'hline', at: 0, color: '--muted', width: 1 },
        { type: 'line', data: zip(d.th, d.a05), color: '--c2', name: 'α = 0.5', fmt: v => fmt(v, 2) },
        { type: 'line', data: zip(d.th, d.a10), color: '--c1', name: 'α = 10', fmt: v => fmt(v, 2) },
        { type: 'points', color: '--c1', data: [[0, v10], [PI, vpi]],
          tip: p => [['--c1', fmt(p[1], 1), p[0] ? 'en θ = π: (α − 1)/√(1 − β)' : 'en θ = 0: 1 − α, el del paper']] },
        { type: 'text', x: 0, y: v10, dx: 10, dy: -9, label: `θ = 0: 1 − α = ${fmt(v10, 0)} (el del paper)` },
        { type: 'text', x: PI, y: vpi, dx: -10, dy: -10, anchor: 'end', label: `θ = π: +${vpi.toFixed(0)}` },
      ],
    };
  }

  // modelo simplificado: log|f| en el plano complejo
  DATA.toy.forEach((d, k) => {
    SPECS['toy' + k] = {
      title: k ? 'Linealizando en θ = π' : 'Modelo del paper, Ec. (4.13)',
      sub: k ? 'm_{eff}^{2} = (α − 1)/√(1 − β), ∂_{a}m* = 0' : 'λ^{2} = 10^{−3}, α = 10, c_{s}^{2} = 0.3, x̄_{0} = 0.6',
      aria: k ? 'Mapa de log|f| linealizando en θ = π: los modos fluidos quedan sobre el eje real y los axiónicos se desplazan.'
              : 'Mapa de log|f| del modelo del paper: diez modos fluidos junto al eje real y cinco axiónicos cada vez más amortiguados.',
      legend: [['heat', '--c1', 'log_{10}|f| (más intenso: más cerca de un cero)'], ['ring', '--c2', 'raíces (Muller)']],
      h: 270, x: { min: d.x[0], max: d.x[1], label: 'Re ω̄' }, y: { min: d.y[0], max: d.y[1], label: 'Im ω̄' },
      layers: [
        { type: 'heat', ...d, vmin: -1.5, vmax: 6 },
        { type: 'points', color: '--c2', hollow: true, r: 4, data: d.roots,
          tip: p => [['--c2', `${fmt(p[0], 2)} ${p[1] < -5e-5 ? '−' : '+'} ${Math.abs(p[1]).toFixed(2)}i`,
            Math.abs(p[1]) < 0.5 && (k === 0 || Math.abs(p[1]) < 1e-6) ? 'modo fluido' : 'modo axiónico']] },
      ],
    };
  });

  // niveles ligados (Tabla IV) y eigenfunciones a 12 kHz
  {
    const d = DATA.levels, L3 = d.levels['3'], L4 = d.levels['4'];
    const all = [...L3.lo, ...L3.hi, ...L4.lo, ...L4.hi, ...L3.table, ...L4.table, 1];
    const lo = Math.min(...all), hi = Math.max(...all), pad = (hi - lo) * 0.12;
    const tl = (n, Ln, c) => [
      { type: 'area', data: d.nu.map((x, i) => [x, Ln.lo[i], Ln.hi[i]]), color: c, opacity: 0.18 },
      { type: 'line', data: zip(d.nu, Ln.mid), color: c, name: `n = ${n}, este cálculo`, fmt: v => v.toFixed(4) },
      { type: 'points', hollow: true, color: c, data: zip(d.nu, Ln.table),
        tip: p => [[c, p[1].toFixed(2), `n = ${n}, Tabla IV a ${p[0]} kHz`], [c, Ln.mid[d.nu.indexOf(p[0])].toFixed(4), 'este cálculo']] },
    ];
    SPECS.levels = {
      title: 'Niveles ligados l = 2 (Tabla IV)', sub: 'la banda es lo que cambia según el interior (plano o Dirichlet, R = 10.5–12 km)',
      aria: 'Cociente −k̃² n²/α² frente a ν_a para n = 3 y 4: las curvas calculadas pasan por los valores de la Tabla IV.',
      legend: [['line', '--c1', 'n = 3'], ['line', '--c2', 'n = 4'], ['area', '--muted', 'banda del interior'], ['ring', '--muted', 'Tabla IV del paper']],
      x: { min: d.nu[0] - 0.3, max: d.nu[d.nu.length - 1] + 0.3, label: 'ν_{a} (kHz)', fmt: v => v + ' kHz' },
      y: { min: lo - pad, max: hi + pad, label: '−k̃_{n}^{2} n^{2} / α^{2}' },
      layers: [{ type: 'hline', at: 1, color: '--muted', dash: '4 3', label: 'hidrogénico, Ec. (23)', lx: d.nu[0] - 0.2 },
        ...tl(3, L3, '--c1'), ...tl(4, L4, '--c2')],
    };
    SPECS.psi = {
      title: 'Eigenfunciones a ν_{a} = 12 kHz', sub: 'como la Fig. 5 del paper; interior «plano», R = 10.5 km',
      aria: 'Potencial ΔV y eigenfunciones de los niveles n = 3 y n = 4 a 12 kHz: el pico de n = 3 está hacia 50 km y n = 4 tiene un nodo hacia 85 km.',
      legend: [['line', '--ink', 'ΔV/m_{a}^{2}, Ec. (22)'], ['line', '--c1', 'ψ, n = 3'], ['line', '--c2', 'ψ, n = 4']],
      x: { min: 0, max: 190, label: 'r (km)', fmt: v => 'r = ' + v.toFixed(0) + ' km' }, y: { min: -0.12, max: 0.06, label: 'energía / m_{a}^{2}' },
      layers: [
        { type: 'vline', at: d.R, color: '--muted', dash: '2 3', label: 'R', ly: -0.112 },
        { type: 'hline', at: d.psi['3'].k2, color: '--c1', width: 0.8 },
        { type: 'hline', at: d.psi['4'].k2, color: '--c2', width: 0.8 },
        { type: 'line', data: zip(d.r, d.dV), color: '--ink', width: 1.6, name: 'ΔV/m_{a}^{2}', fmt: v => fmt(v, 4) },
        { type: 'line', data: zip(d.r, d.psi['3'].y), color: '--c1' },
        { type: 'line', data: zip(d.r, d.psi['4'].y), color: '--c2' },
      ],
    };
  }

  // resonancia del modo f y las 52 frecuencias
  {
    const d = DATA.modes;
    const all = [...d.full, ...d.zl].filter(v => v > 0);
    const ymin = 10 ** Math.floor(Math.log10(Math.min(...all)));
    const yt = [];
    for (let e = Math.log10(ymin); e <= 0; e += 2) yt.push([10 ** e, sci(10 ** e)]);
    SPECS.res = {
      title: 'Modo f, l = 2: la resonancia', sub: 'politropo de 1.4 M_{⊙}; el mínimo de la onda saliente marca el modo',
      aria: `Amplitud de la onda saliente frente a la frecuencia: mínimo en ${d.lmi}, el valor de Lindblom, Mendell e Ipser; Cowling da ${d.cow.toFixed(3)}; con la errata de signo, ${d.zlmin.toFixed(3)}.`,
      legend: [['line', '--c1', 'métrica dinámica (sistema extendido, sin axión)'], ['dash', '--c3', 'ídem, con el signo impreso por Zhao y Lattimer'], ['line', '--c2', 'Cowling (paso 5)']],
      x: { min: 0.6, max: 2.0, label: 'ω / √(M/R^{3})', fmt: v => 'ω = ' + v.toFixed(4) },
      y: { min: ymin, max: 2, log: true, ticks: yt, label: '|A_{out}|^{2} normalizada' }, ml: 50,
      layers: [
        { type: 'vline', at: d.lmi, color: '--muted', dash: '2 3', label: `LMI (1997): ${d.lmi}`, ly: ymin * 30 },
        { type: 'vline', at: d.cow, color: '--c2', width: 2, label: `Cowling: ${d.cow.toFixed(3)}`, ly: 0.12 },
        { type: 'line', data: zip(d.scan, d.zl), color: '--c3', dash: '5 3', name: 'con la errata', fmt: v => v.toExponential(1) },
        { type: 'line', data: zip(d.scan, d.full), color: '--c1', name: 'métrica dinámica', fmt: v => v.toExponential(1) },
        { type: 'text', x: d.zlmin, y: ymin * 4, dx: -6, anchor: 'end', label: `errata: ${d.zlmin.toFixed(3)}` },
      ],
    };
    const cols = ['--c2', '--c1', '--c1'], names = ['Cowling vs. YK', 'dinámica vs. YK', 'dinámica vs. LMI'];
    let x0 = 0;
    const layers = [{ type: 'hline', at: 0, color: '--muted', width: 1 }], ticks = [];
    let ext = 0;
    d.groups.forEach((gr, k) => {
      const data = gr.pts.map((p, i) => [x0 + i, p[0], p[1], gr.tags[i], p[2], p[3]]);
      ext = Math.max(ext, ...gr.pts.map(p => Math.abs(p[0]) + p[1]));
      layers.push({ type: 'points', err: true, hollow: k === 2, color: cols[k], data, r: 3.5,
        tip: p => [[cols[k], (p[1] >= 0 ? '+' : '−') + Math.abs(p[1]).toFixed(3) + ' %', names[k]], [null, p[3], ''],
          [null, `${p[4].toPrecision(5)} frente a ${p[5]}`, '']] });
      ticks.push([x0 + (gr.pts.length - 1) / 2, names[k]]);
      x0 += gr.pts.length + 4;
    });
    SPECS.cmp = {
      title: '52 frecuencias publicadas', sub: 'diferencia relativa; la barra es la precisión con que se publicó cada valor',
      aria: 'Diferencias relativas con 52 frecuencias publicadas: casi todas dentro de su precisión; las de Yoshida y Kojima en relatividad completa, dentro del 0.5 %.',
      x: { min: -2, max: x0 - 2, ticks, nogrid: true, label: '' }, mb: 30,
      y: { min: -ext * 1.15, max: ext * 1.15, label: 'este cálculo − referencia (%)' }, layers,
    };
  }

  // fondo con axión
  {
    const d = DATA.background, cols = ['--c2', '--c1'];
    const layers = [{ type: 'vline', at: 1, color: '--muted', dash: '2 3', label: 'superficie', ly: 0.97 }];
    d.bg.forEach((b, k) => layers.push({ type: 'area', data: d.x.map((x, i) => [x, 0, b.carrier[i]]), color: cols[k], opacity: 0.16 }));
    d.bg.forEach((b, k) => {
      layers.push({ type: 'line', data: zip(d.x, b.theta), color: cols[k], name: `θ_{0}/π, ${b.nu} kHz`, fmt: v => v.toFixed(3) });
      layers.push({ type: 'line', data: zip(d.x, b.carrier), color: cols[k], width: 0.01, name: `acoplo, ${b.nu} kHz`, fmt: v => v.toFixed(2) });
    });
    SPECS.bg = {
      title: 'Fondo con axión: dónde actúan los términos nuevos', sub: 'esquemático: politropo n = 1 de 1.4 M_{⊙} y 12 km, sin retroacción sobre la métrica',
      aria: 'θ₀/π frente a r/R para 5 y 15 kHz, con sombras donde el acoplo es grande: a 15 kHz, una capa entre 0.65 R y la superficie; a 5 kHz, media estrella.',
      legend: d.bg.map((b, k) => ['line', cols[k], `θ_{0}/π, ν_{a} = ${b.nu} kHz`]).concat([['area', '--muted', 'n_{b}|∂m*/∂θ| normalizado']]),
      x: { min: 0, max: 1.6, label: 'r / R', fmt: v => 'r = ' + v.toFixed(2) + ' R' }, y: { min: 0, max: 1.05, label: 'θ_{0}/π  ·  acoplo (sombra)' }, layers,
    };
  }

  // con axión en Cowling (paso 9): desplazamientos de las frecuencias y amortiguamiento
  if (DATA.axcow) {
    const d = DATA.axcow, lab = { f: 'f', p1: 'p_{1}', p2: 'p_{2}' };
    const X = { min: 2.3, max: 15.5, label: 'ν_{a} (kHz)', fmt: v => 'ν_a = ' + v + ' kHz', name: '' };
    d.modes.forEach(md => {
      const all = [...md.bg, ...md.cp, ...md.pbg.map(p => p[1]), ...md.pcp.map(p => p[1]), 0];
      const lo = Math.min(...all), hi = Math.max(...all), pad = (hi - lo) * 0.1;
      const pct = v => (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(2) + ' %';
      SPECS['ax_' + md.name] = {
        title: `Modo ${lab[md.name]}`, sub: 'desplazamiento respecto a la estrella sin axión', h: 250,
        aria: `Desplazamiento del modo ${md.name} frente a la masa del axión: este cálculo y la secuela coinciden dentro de medio punto.`,
        legend: [['line', '--c2', 'solo fondo (Tabla I)'], ['line', '--c1', 'con δθ (Tabla II)'], ['ring', '--muted', 'secuela']],
        x: X, y: { min: lo - pad, max: hi + pad, label: 'ν/ν_{0} − 1 (%)' },
        layers: [
          { type: 'hline', at: 0, color: '--muted', width: 1 },
          { type: 'line', data: zip(d.nus, md.bg), color: '--c2', name: 'solo fondo', fmt: pct },
          { type: 'line', data: zip(d.nus, md.cp), color: '--c1', name: 'con δθ', fmt: pct },
          { type: 'points', hollow: true, color: '--c2', r: 3.5, data: md.pbg, tip: p => [['--c2', pct(p[1]), `secuela, solo fondo, ${p[0]} kHz`]] },
          { type: 'points', hollow: true, color: '--c1', r: 3.5, data: md.pcp, tip: p => [['--c1', pct(p[1]), `secuela, con δθ, ${p[0]} kHz`]] },
        ],
      };
    });
    const cols = { p1: '--c3', p2: '--ink' };
    const layers = [];
    d.modes.filter(md => md.name !== 'f').forEach(md => {
      const ours = d.nus.map((nu, i) => [nu, md.im[i]]).filter(p => p[1] > 1e-3);
      layers.unshift({ type: 'points', color: cols[md.name], hollow: true, r: 6, data: md.pim,
        tip: p => [[cols[md.name], p[1].toPrecision(2) + ' Hz', `${md.name}, secuela, ${p[0]} kHz`]] });
      layers.push({ type: 'points', color: cols[md.name], r: 3.5, data: ours,
        tip: p => [[cols[md.name], p[1].toPrecision(2) + ' Hz', `${md.name}, este cálculo, ${p[0]} kHz`],
          [null, (1e3 / (2 * PI * p[1])).toFixed(0) + ' ms', 'tiempo de amortiguamiento']] });
    });
    SPECS.ax_im = {
      title: 'Amortiguamiento por emisión de axiones', sub: '|Im ν|; τ_{a} = 1/(2π|Im ν|), de 12 ms a varios segundos', h: 250,
      aria: 'Amortiguamiento de p1 y p2 frente a la masa del axión: mismo orden de magnitud que la secuela, con el máximo de unos 13 Hz entre 5 y 6 kHz.',
      legend: [['dot', '--c3', 'p_{1}'], ['dot', '--ink', 'p_{2}'], ['dot', '--muted', 'este cálculo'], ['ring', '--muted', 'secuela']],
      x: X, y: { min: 0.03, max: 30, log: true, ticks: [[0.1, '0.1'], [1, '1'], [10, '10']], label: '|Im ν| (Hz)' }, layers,
    };
  }

  // etapa B (paso 10): con axión y métrica dinámica, frente a Cowling
  if (DATA.axgr && DATA.axcow) {
    const d = DATA.axgr, c = DATA.axcow, lab = { f: 'f', p1: 'p_{1}', p2: 'p_{2}' };
    const X = { min: 2.3, max: 15.5, label: 'ν_{a} (kHz)', fmt: v => 'ν_a = ' + v + ' kHz', name: '' };
    const pct = v => (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(2) + ' %';
    d.modes.forEach((md, i) => {
      const cw = c.modes[i];
      const all = [...md.cp, ...md.bg, ...cw.cp, 0];
      const lo = Math.min(...all), hi = Math.max(...all), pad = (hi - lo) * 0.1;
      SPECS['gr_' + md.name] = {
        title: `Modo ${lab[md.name]}`, sub: `sin axión: ${d.gr0[i][0].toFixed(3)} kHz con la métrica dinámica, ${d.cow0[i].toFixed(3)} en Cowling`, h: 250,
        aria: `Desplazamiento del modo ${md.name} con la métrica dinámica frente a la masa del axión, comparado con Cowling.`,
        legend: [['line', '--c1', 'métrica dinámica, con δθ (nuevo)'], ['line', '--c2', 'métrica dinámica, solo fondo'], ['dash', '--muted', 'Cowling, con δθ (paso 9)']],
        x: X, y: { min: lo - pad, max: hi + pad, label: 'ν/ν_{0} − 1 (%)' },
        layers: [
          { type: 'hline', at: 0, color: '--muted', width: 1 },
          { type: 'line', data: zip(c.nus, cw.cp), color: '--muted', dash: '5 3', name: 'Cowling', fmt: pct },
          { type: 'line', data: zip(d.nus, md.bg), color: '--c2', name: 'solo fondo', fmt: pct },
          { type: 'line', data: zip(d.nus, md.cp), color: '--c1', name: 'con δθ', fmt: pct },
          { type: 'points', color: '--c1', r: 3.5, data: zip(d.nus, md.cp),
            tip: p => [['--c1', pct(p[1]), `con δθ, ${p[0]} kHz`], [null, md.re[d.nus.indexOf(p[0])].toFixed(3) + ' kHz', '']] },
        ],
      };
    });
    const cols = { f: '--c1', p1: '--c3', p2: '--ink' }, layers = [];
    c.modes.filter(md => md.name !== 'f').forEach(md => {
      const pts = c.nus.map((nu, i) => [nu, md.im[i]]).filter(p => p[1] > 1e-3);
      layers.push({ type: 'points', color: cols[md.name], hollow: true, r: 6, data: pts,
        tip: p => [[cols[md.name], p[1].toPrecision(2) + ' Hz', `${md.name}, Cowling (solo axiones), ${p[0]} kHz`]] });
    });
    d.modes.forEach(md => {
      const pts = d.nus.map((nu, i) => [nu, md.im[i]]).filter(p => p[1] > 1e-3);
      layers.push({ type: 'points', color: cols[md.name], r: 3.5, data: pts,
        tip: p => [[cols[md.name], p[1].toPrecision(3) + ' Hz', `${md.name}, métrica dinámica, ${p[0]} kHz`],
          [null, (1e3 / (2 * PI * p[1])).toPrecision(3) + ' ms', 'tiempo de amortiguamiento']] });
    });
    SPECS.gr_im = {
      title: 'Amortiguamiento total', sub: 'ondas gravitacionales y axiones a la vez; anillos: solo axiones en Cowling', h: 270,
      aria: 'Amortiguamiento de f, p1 y p2 con la métrica dinámica frente a la masa del axión, comparado con el de Cowling.',
      legend: [['dot', '--c1', 'f'], ['dot', '--c3', 'p_{1}'], ['dot', '--ink', 'p_{2}'], ['dot', '--muted', 'métrica dinámica'], ['ring', '--muted', 'Cowling']],
      x: X, y: { min: 1e-3, max: 30, log: true, ticks: [[0.001, '0.001'], [0.01, '0.01'], [0.1, '0.1'], [1, '1'], [10, '10']], label: '|Im ν| (Hz)' }, layers,
    };
  }

  // ------------------------------------------------------------------ montaje y redibujado
  const boxes = [...document.querySelectorAll('[data-chart]')];
  const widths = new Map();
  const draw = b => { widths.set(b, b.clientWidth); render(b, SPECS[b.dataset.chart]); };
  boxes.forEach(draw);
  const ro = new ResizeObserver(es => es.forEach(e => {
    if (Math.abs(e.target.clientWidth - (widths.get(e.target) || 0)) > 2) draw(e.target);
  }));
  boxes.forEach(b => ro.observe(b));
  const redrawAll = () => boxes.forEach(draw);
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', redrawAll);
  new MutationObserver(redrawAll).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

  // tabla de las 52 frecuencias
  const tb = document.getElementById('freqtable');
  if (tb) {
    const names = ['Cowling vs. Yoshida–Kojima', 'Métrica dinámica vs. Yoshida–Kojima', 'Métrica dinámica vs. Lindblom–Mendell–Ipser'];
    DATA.modes.groups.forEach((gr, k) => gr.pts.forEach((p, i) => {
      const tr = htmlEl('tr', '', tb);
      for (const [v, c] of [[i ? '' : names[k], ''], [gr.tags[i], ''], [p[2].toPrecision(5), 'num'], [String(p[3]), 'num'],
        [(p[0] >= 0 ? '+' : '−') + Math.abs(p[0]).toFixed(3), 'num'], ['±' + p[1].toFixed(3), 'num']])
        htmlEl('td', c, tr).textContent = v;
    }));
  }
})();
