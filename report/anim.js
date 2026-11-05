// Animaciones estilo 3Blue1Brown / Manim:
// 1) "De la acción a los modos, en seis ecuaciones" (escenas vectoriales en vivo + ecuaciones).
// 2) "Anatomía de un modo" (datos reales de los pasos 9 y 10 con retícula de espacio-tiempo, vectores de velocidad,
//    capa de transición, nodos radiales y osciloscopio de ringdown en tiempo real).
// Sin bibliotecas externas: Canvas 2D y CSS. Con prefers-reduced-motion, la animación empieza en pausa.
(() => {
  const DATA = JSON.parse(document.getElementById('figdata').textContent);
  const REDUCE = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const PAL = {
    bg: [14, 18, 28], star: [28, 43, 74], cond: [247, 217, 111],
    gwp: [88, 196, 221], gwn: [154, 114, 172],
    axp: [247, 217, 111], axn: [252, 98, 85],
    blue: '#58c4dd', yellow: '#f7d96f', green: '#83c167', red: '#fc6255', purple: '#9a72ac',
    ink: '#e6edf6', muted: '#8b97ab', grid: 'rgba(88, 196, 221, 0.14)'
  };

  function arrow(ctx, x1, y1, x2, y2, color, size = 6, lw = 1.6) {
    const dx = x2 - x1, dy = y2 - y1, len = Math.hypot(dx, dy);
    if (len < 1e-3) return;
    const ux = dx / len, uy = dy / len;
    ctx.save();
    ctx.strokeStyle = color;
    ctx.fillStyle = color;
    ctx.lineWidth = lw;
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2 - ux * size * 0.7, y2 - uy * size * 0.7);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x2, y2);
    ctx.lineTo(x2 - ux * size - uy * size * 0.52, y2 - uy * size + ux * size * 0.52);
    ctx.lineTo(x2 - ux * size + uy * size * 0.52, y2 - uy * size - ux * size * 0.52);
    ctx.closePath();
    ctx.fill();
    ctx.restore();
  }

  function glowStroke(ctx, color, width, glow = 8) {
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.shadowColor = color;
    ctx.shadowBlur = glow;
    ctx.stroke();
    ctx.shadowBlur = 0;
  }

  // ------------------------------------------------------------------ 1. anatomía de un modo
  function modeAnim(root) {
    const gr = DATA.axgr, cw = DATA.axcow;
    if (!gr || !cw) return;
    const cv = root.querySelector('canvas'), ctx = cv.getContext('2d');
    const off = document.createElement('canvas'), octx = off.getContext('2d');
    const ui = {
      modes: root.querySelectorAll('[data-mode]'), metric: root.querySelectorAll('[data-metric]'),
      slider: root.querySelector('input[type=range]'), nu: root.querySelector('[data-out=nua]'),
      freq: root.querySelector('[data-out=freq]'), damp: root.querySelector('[data-out=damp]'),
      rad: root.querySelector('[data-out=rad]'), play: root.querySelector('[data-play]'),
    };
    const st = { mode: 2, idx: 1 + gr.nus.indexOf(5), dyn: true, playing: !REDUCE, t: 0, visible: true };
    ui.slider.max = gr.nus.length;
    ui.slider.value = st.idx;

    function current() {
      const i = st.mode, j = st.idx - 1, nua = j >= 0 ? gr.nus[j] : null;
      let re, im;
      if (st.dyn) {
        re = j >= 0 ? gr.modes[i].re[j] : gr.gr0[i][0];
        im = j >= 0 ? gr.modes[i].im[j] : gr.gr0[i][1];
      } else {
        re = j >= 0 ? gr.cow0[i] * (1 + cw.modes[i].cp[j] / 100) : gr.cow0[i];
        im = j >= 0 ? cw.modes[i].im[j] : 0;
      }
      const prof = j >= 0 ? cw.prof[j] : null;
      return { i, nua, re, im, prof, radiates: nua !== null && re > nua };
    }

    let W, H, R, cx, cy, S = 3, grid, img;
    function layout() {
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      W = Math.max(300, root.querySelector('.stage-canvas').clientWidth);
      H = Math.round(Math.min(450, Math.max(310, W * 0.53)));
      cv.width = W * dpr; cv.height = H * dpr; cv.style.height = H + 'px';
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      R = Math.min(H * 0.235, W * 0.13); cx = W / 2; cy = H / 2 - 6;
      off.width = Math.ceil(W / S); off.height = Math.ceil(H / S);
      const n = off.width * off.height;
      grid = { rho: new Float32Array(n), c2: new Float32Array(n), side: new Int8Array(n) };
      for (let y = 0; y < off.height; y++) for (let x = 0; x < off.width; x++) {
        const dx = x * S - cx, dy = y * S - cy, k = y * off.width + x;
        const rho = Math.hypot(dx, dy);
        grid.rho[k] = rho; grid.c2[k] = rho > 0 ? (dy * dy) / (rho * rho) : 1;
        grid.side[k] = dx < 0 ? -1 : 1;
      }
      img = octx.createImageData(off.width, off.height);
    }

    // perfil radial del desplazamiento (estilizado): f sin nodos, p1 uno (x≈0.67), p2 dos (x≈0.46, 0.89)
    const xi = [x => x, x => x * (1 - 2.2 * x * x) / -1.2, x => x * (1 - 6 * x * x + 6 * x ** 4)];
    const nodes = [[], [0.674], [0.459, 0.888]];
    const mix = (a, b, t) => a + (b - a) * t;

    function transitionRadius(prof) {
      if (!prof) return null;
      for (let i = 0; i < prof.length - 1; i++) {
        if (prof[i] >= 0.5 && prof[i + 1] < 0.5) {
          const u = (prof[i] - 0.5) / (prof[i] - prof[i + 1]);
          return (i + u) * 0.05;
        }
      }
      return 0.72;
    }

    function draw() {
      const m = current();
      const nuRef = 1.961;
      const Lref = 2.4 * R;
      const lamGW = Lref * nuRef / m.re;
      let lamA = null, ell = null;
      if (m.nua !== null) {
        if (m.radiates) lamA = Lref * nuRef / Math.sqrt(m.re ** 2 - m.nua ** 2);
        else ell = Lref * nuRef / (2 * Math.PI * Math.sqrt(m.nua ** 2 - m.re ** 2));
      }
      const ph = st.t * 2 * Math.PI / 1.4;
      const amp = 0.075, cosw = Math.cos(ph), sinw = Math.sin(ph);
      const d = img.data, prof = m.prof;

      for (let k = 0; k < grid.rho.length; k++) {
        const rho = grid.rho[k], c2 = grid.c2[k], P2 = 1.5 * c2 - 0.5;
        const rb = R * (1 + amp * xi[m.i](1) * P2 * cosw);
        let c = PAL.bg, a = 0, col = null;
        const o = k * 4;
        if (rho <= rb) {
          const x = rho / R;
          let th = 0;
          if (prof) {
            const u = Math.min(x / 0.05, prof.length - 1.001), i0 = Math.floor(u);
            th = mix(prof[i0], prof[i0 + 1], u - i0);
          }
          const tone = 0.15 + 0.85 * th;
          const disp = 0.5 + 0.5 * xi[m.i](Math.min(x, 1)) * P2 * cosw;
          for (let q = 0; q < 3; q++) {
            d[o + q] = mix(mix(PAL.star[q], PAL.cond[q], th ** 1.5 * 0.78), 255, 0.12 * disp) * (0.7 + 0.3 * tone);
          }
          d[o + 3] = 255;
          continue;
        }
        const s = rho - R;
        if (grid.side[k] < 0) {
          if (st.dyn) {
            const v = (1 - c2) * Math.cos(ph - 2 * Math.PI * s / lamGW) * Math.min(1, 1.6 * R / rho);
            col = v > 0 ? PAL.gwp : PAL.gwn; a = Math.min(1, Math.abs(v)) * 0.85;
          }
        } else if (m.nua !== null) {
          let v = 0;
          if (lamA) v = P2 * Math.cos(ph - 2 * Math.PI * s / lamA) * Math.min(1, 1.6 * R / rho);
          else v = P2 * cosw * Math.exp(-s / ell);
          if (prof) {
            const x = rho / R, u = Math.min(x / 0.05, prof.length - 1.001), i0 = Math.floor(u);
            const th = mix(prof[i0], prof[i0 + 1], u - i0);
            c = PAL.bg.map((v0, q) => mix(v0, PAL.cond[q], 0.35 * th));
          }
          col = v > 0 ? PAL.axp : PAL.axn; a = Math.min(1, Math.abs(v)) * 0.82;
        }
        for (let q = 0; q < 3; q++) d[o + q] = col ? mix(c[q], col[q], a) : c[q];
        d[o + 3] = 255;
      }
      octx.putImageData(img, 0, 0);
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(off, 0, 0, off.width * S, off.height * S);

      // Retícula de espacio-tiempo 3Blue1Brown (en la izquierda se curva con h+ si métrica dinámica, y queda rígida en Cowling)
      ctx.save();
      ctx.strokeStyle = st.dyn ? 'rgba(88, 196, 221, 0.19)' : 'rgba(139, 151, 171, 0.12)';
      ctx.lineWidth = 1;
      for (let gx = -5; gx <= -1; gx++) {
        const x0 = cx + gx * R * 0.62;
        if (x0 < 12) continue;
        ctx.beginPath();
        for (let gy = -24; gy <= 24; gy++) {
          const y0 = cy + gy * (H / 52);
          const dx = x0 - cx, dy = y0 - cy, rho = Math.hypot(dx, dy);
          let wx = x0, wy = y0;
          if (st.dyn && rho > R * 1.05) {
            const c2 = (dy * dy) / (rho * rho);
            const h = 4.5 * (1 - c2) * Math.cos(ph - 2 * Math.PI * (rho - R) / lamGW) * Math.min(1, 1.5 * R / rho);
            wx += h * (dx / rho);
            wy -= h * (dy / rho);
          }
          gy === -24 ? ctx.moveTo(wx, wy) : ctx.lineTo(wx, wy);
        }
        ctx.stroke();
      }
      for (let gy = -3; gy <= 3; gy++) {
        const y0 = cy + gy * R * 0.58;
        ctx.beginPath();
        for (let gx = -32; gx <= 0; gx++) {
          const x0 = cx + gx * (W / 68);
          const dx = x0 - cx, dy = y0 - cy, rho = Math.hypot(dx, dy);
          if (rho < R * 1.04) { ctx.moveTo(x0, y0); continue; }
          let wx = x0, wy = y0;
          if (st.dyn) {
            const c2 = (dy * dy) / (rho * rho);
            const h = 4.5 * (1 - c2) * Math.cos(ph - 2 * Math.PI * (rho - R) / lamGW) * Math.min(1, 1.5 * R / rho);
            wx += h * (dx / rho);
            wy -= h * (dy / rho);
          }
          gx === -32 ? ctx.moveTo(wx, wy) : ctx.lineTo(wx, wy);
        }
        ctx.stroke();
      }
      ctx.restore();

      // Divisor central vertical sutil
      ctx.save();
      ctx.setLineDash([3, 5]);
      ctx.strokeStyle = 'rgba(230, 237, 246, 0.18)';
      ctx.beginPath();
      ctx.moveTo(cx, 12); ctx.lineTo(cx, cy - R * 1.15);
      ctx.moveTo(cx, cy + R * 1.15); ctx.lineTo(cx, H - 12);
      ctx.stroke();
      ctx.restore();

      // Capas del fluido + superficie con brillo Manim
      for (const r0 of [0.25, 0.5, 0.75, 1]) {
        ctx.beginPath();
        for (let a2 = 0; a2 <= 72; a2++) {
          const t = a2 / 72 * 2 * Math.PI, c = Math.cos(t), P2 = 1.5 * c * c - 0.5;
          const rr = R * (r0 + amp * xi[m.i](r0) * P2 * cosw);
          const X = cx + rr * Math.sin(t), Y = cy - rr * c;
          a2 ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y);
        }
        if (r0 === 1) {
          glowStroke(ctx, 'rgba(230,237,246,0.95)', 2.2, 6);
        } else {
          ctx.strokeStyle = 'rgba(14,18,28,0.45)';
          ctx.lineWidth = 1;
          ctx.stroke();
        }
      }

      // Nodos radiales de p1 / p2 en cian punteado (#58c4dd)
      ctx.save();
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = PAL.blue;
      ctx.lineWidth = 1.4;
      for (const rn of nodes[m.i]) {
        ctx.beginPath();
        ctx.arc(cx, cy, R * rn, 0, 2 * Math.PI);
        ctx.stroke();
      }
      ctx.restore();

      // Vectores de velocidad cuadrupolar del fluido dentro de la estrella (estilo 3B1B)
      const vphase = -sinw;
      if (Math.abs(vphase) > 0.12) {
        const rSample = m.i === 0 ? 0.62 : m.i === 1 ? 0.84 : 0.68;
        const sgn = Math.sign(xi[m.i](rSample)) * vphase;
        for (const ang of [0, Math.PI / 2, Math.PI, 1.5 * Math.PI]) {
          const c = Math.cos(ang), s = Math.sin(ang), P2 = 1.5 * c * c - 0.5;
          const rBase = R * rSample;
          const x1 = cx + rBase * s, y1 = cy - rBase * c;
          const dr = 11 * P2 * sgn;
          arrow(ctx, x1, y1, x1 + dr * s, y1 - dr * c, 'rgba(230,237,246,0.85)', 4.5, 1.3);
        }
      }

      // Anillo y llamada vectorial a la capa de transición (θ₀ ≈ π/2)
      const xTrans = transitionRadius(prof);
      if (xTrans !== null && W > 520) {
        const rT = R * xTrans;
        ctx.save();
        ctx.setLineDash([2, 3]);
        ctx.strokeStyle = PAL.yellow;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.arc(cx, cy, rT, -0.7, 0.7);
        ctx.stroke();
        ctx.restore();

        const tx = cx + rT * Math.cos(-0.42), ty = cy + rT * Math.sin(-0.42);
        const lx = cx + R * 1.48, ly = cy - R * 0.92;
        arrow(ctx, lx - 4, ly + 4, tx + 3, ty - 2, PAL.yellow, 5.5, 1.4);
        ctx.font = '500 11.5px "IBM Plex Mono", monospace';
        ctx.fillStyle = PAL.yellow;
        ctx.textAlign = 'left';
        ctx.fillText('capa de transición (∂m*/∂θ ≠ 0)', lx, ly);
      }

      // Rótulos superiores estilo 3B1B
      ctx.font = '600 13px "IBM Plex Sans", system-ui, sans-serif';
      ctx.fillStyle = PAL.blue; ctx.textAlign = 'left';
      ctx.fillText('← Ondas gravitacionales (h₊)', 14, 24);
      ctx.fillStyle = PAL.yellow; ctx.textAlign = 'right';
      ctx.fillText('Radiación axiónica (δθ) →', W - 14, 24);

      ctx.font = '400 12px "IBM Plex Sans", system-ui, sans-serif';
      ctx.fillStyle = PAL.muted;
      ctx.textAlign = 'left';
      ctx.fillText(st.dyn ? 'Métrica dinámica: retícula oscila (∝ sin²ϑ)' : 'Cowling: métrica congelada (retícula fija)', 14, 42);
      ctx.textAlign = 'right';
      ctx.fillText(m.nua === null ? 'Sin condensado axiónico' : m.radiates ? `ν (${m.re.toFixed(2)}) > νₐ (${m.nua}): onda escapa` : `ν (${m.re.toFixed(2)}) < νₐ (${m.nua}): nube ligada`, W - 14, 42);

      // Etiqueta de nodos radiales debajo de la estrella
      const nodeLabel = ['modo f: 0 nodos radiales', 'modo p₁: 1 nodo radial (círculo azul)', 'modo p₂: 2 nodos radiales'][m.i];
      ctx.textAlign = 'center';
      ctx.font = '500 11.5px "IBM Plex Mono", monospace';
      ctx.fillStyle = nodes[m.i].length ? PAL.blue : PAL.ink;
      ctx.fillText(nodeLabel, cx, cy + R * 1.28);

      // Mini-osciloscopio Manim en esquina inferior izquierda: amortiguamiento temporal e^(-t/τ) cos(ωt)
      if (W > 560) {
        const bx = 14, by = H - 68, bw = 175, bh = 44;
        ctx.save();
        ctx.fillStyle = 'rgba(20, 27, 39, 0.86)';
        ctx.strokeStyle = '#2c3850';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.roundRect ? ctx.roundRect(bx, by, bw, bh, 5) : ctx.rect(bx, by, bw, bh);
        ctx.fill(); ctx.stroke();

        // tasa visual de caída (escalada logarítmicamente a partir de m.im real)
        const decay = m.im > 1e-4 ? Math.min(2.6, 0.32 * Math.log10(1 + m.im * 18)) : 0;
        const midY = by + bh / 2, ampY = bh * 0.36;

        // envolvente ±e^(-t/τ)
        if (decay > 0.01) {
          ctx.setLineDash([2, 2]);
          ctx.strokeStyle = 'rgba(252, 98, 85, 0.55)';
          for (const sgn of [-1, 1]) {
            ctx.beginPath();
            for (let u = 0; u <= bw - 16; u += 3) {
              const env = Math.exp(-decay * (u / (bw - 16)));
              const px = bx + 8 + u, py = midY + sgn * ampY * env;
              u === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
            }
            ctx.stroke();
          }
          ctx.setLineDash([]);
        }

        // traza oscilante
        ctx.beginPath();
        for (let u = 0; u <= bw - 16; u += 1.5) {
          const frac = u / (bw - 16);
          const env = Math.exp(-decay * frac);
          const y = midY - ampY * env * Math.cos(frac * 6.5 * Math.PI - ph);
          u === 0 ? ctx.moveTo(bx + 8 + u, y) : ctx.lineTo(bx + 8 + u, y);
        }
        glowStroke(ctx, m.radiates ? PAL.yellow : PAL.blue, 1.6, 4);

        ctx.font = '500 10px "IBM Plex Mono", monospace';
        ctx.fillStyle = PAL.muted;
        ctx.textAlign = 'left';
        ctx.fillText(decay > 0.01 ? 'ringdown: e^(-t/τ) cos(ωt)' : 'sin amortiguamiento (τ = ∞)', bx + 6, by - 5);
        ctx.restore();
      }

      if (W > 720) {
        ctx.font = '400 11.5px "IBM Plex Sans", system-ui, sans-serif';
        ctx.fillStyle = PAL.muted;
        ctx.textAlign = 'right';
        ctx.fillText('distancias exteriores comprimidas · proporción de longitudes de onda real', W - 14, H - 12);
      }
    }

    function readout() {
      const m = current(), lab = ['f', 'p₁', 'p₂'][m.i];
      ui.nu.textContent = m.nua === null ? 'sin axión' : `νₐ = ${m.nua} kHz`;
      ui.freq.textContent = `${lab}: ${m.re.toFixed(3)} kHz`;
      if (m.im > 1e-3) {
        const tau = 1 / (2 * Math.PI * m.im);
        ui.damp.textContent = `|Im ν| = ${m.im < 0.1 ? m.im.toPrecision(2) : m.im.toFixed(2)} Hz · τ = ${tau >= 1 ? tau.toFixed(2) + ' s' : (tau * 1e3).toPrecision(3) + ' ms'}`;
      } else ui.damp.textContent = 'sin amortiguamiento';
      ui.rad.textContent = m.nua === null ? '' : m.radiates ? 'radia axiones' : 'no radia axiones';
      ui.modes.forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.mode === st.mode)));
      ui.metric.forEach(b => b.setAttribute('aria-pressed', String((b.dataset.metric === 'dyn') === st.dyn)));
      ui.play.textContent = st.playing ? 'Pausa' : 'Animar';
      ui.play.setAttribute('aria-pressed', String(st.playing));
    }

    let last = null;
    function frame(ts) {
      if (last !== null && st.playing) st.t += Math.min(0.05, (ts - last) / 1000);
      last = ts;
      if (st.visible) draw();
      if (st.playing && st.visible) requestAnimationFrame(frame); else last = null;
    }
    const kick = () => { readout(); if (!st.playing || !st.visible) draw(); else if (last === null) requestAnimationFrame(frame); };
    ui.modes.forEach(b => b.addEventListener('click', () => { st.mode = +b.dataset.mode; kick(); }));
    ui.metric.forEach(b => b.addEventListener('click', () => { st.dyn = b.dataset.metric === 'dyn'; kick(); }));
    ui.slider.addEventListener('input', () => { st.idx = +ui.slider.value; kick(); });
    ui.play.addEventListener('click', () => { st.playing = !st.playing; readout(); if (st.playing) requestAnimationFrame(frame); });
    new IntersectionObserver(es => { st.visible = es[0].isIntersecting; if (st.visible && st.playing && last === null) requestAnimationFrame(frame); })
      .observe(cv);
    let lastW = 0;
    new ResizeObserver(() => {
      const w = root.querySelector('.stage-canvas').clientWidth;
      if (Math.abs(w - lastW) > 1) { lastW = w; layout(); draw(); }
    }).observe(root.querySelector('.stage-canvas'));
    layout(); readout(); draw();
    if (st.playing) requestAnimationFrame(frame);
  }

  // ------------------------------------------------------------------ 2. cadena de ecuaciones (6 escenas Manim)
  function stepper(root) {
    const steps = [...root.querySelectorAll('.step')], dots = [...root.querySelectorAll('.rail button:not([data-chain-play])')];
    const prev = root.querySelector('[data-prev]'), next = root.querySelector('[data-next]');
    const playBtn = root.querySelector('[data-chain-play]');
    const cv = root.querySelector('canvas.chain-cv');
    const ctx = cv ? cv.getContext('2d') : null;
    let k = 0, t = 0, stepAnim = 1, playing = !REDUCE, visible = true, last = null, W = 600, H = 210;

    function layout() {
      if (!cv) return;
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      W = Math.max(280, cv.parentElement.clientWidth);
      H = W < 520 ? 195 : 215;
      cv.width = W * dpr; cv.height = H * dpr; cv.style.height = H + 'px';
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    function drawBackgroundGrid() {
      ctx.fillStyle = '#0b0f18';
      ctx.fillRect(0, 0, W, H);
      ctx.strokeStyle = 'rgba(88, 196, 221, 0.06)';
      ctx.lineWidth = 1;
      const step = 28;
      for (let x = (W % step) / 2; x < W; x += step) {
        ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, H); ctx.stroke();
      }
      for (let y = (H % step) / 2; y < H; y += step) {
        ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(W, y); ctx.stroke();
      }
    }

    // Escena 0: Acción (Gravedad + Fluido + Axión conectados por m*(a))
    function scene0(ph) {
      const y = H * 0.52, x1 = W * 0.2, x2 = W * 0.5, x3 = W * 0.8;
      // 1. Gravedad: malla curvada en azul Manim
      for (let i = -2; i <= 2; i++) {
        ctx.beginPath();
        for (let u = -38; u <= 38; u += 4) {
          const r2 = u * u + (i * 13) ** 2;
          const dip = 14 * Math.exp(-r2 / 900) * (1 + 0.15 * Math.sin(ph));
          const px = x1 + u, py = y + i * 12 + dip;
          u === -38 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
        }
        glowStroke(ctx, PAL.blue, 1.5, 4);
      }
      // 2. Fluido: núcleo bariónico verde Manim
      const rCore = 28 + 2.5 * Math.sin(ph * 1.3);
      const grad = ctx.createRadialGradient(x2, y, 4, x2, y, rCore * 1.25);
      grad.addColorStop(0, 'rgba(131, 193, 103, 0.85)');
      grad.addColorStop(1, 'rgba(131, 193, 103, 0.05)');
      ctx.fillStyle = grad;
      ctx.beginPath(); ctx.arc(x2, y, rCore * 1.25, 0, 2 * Math.PI); ctx.fill();
      ctx.beginPath(); ctx.arc(x2, y, rCore, 0, 2 * Math.PI);
      glowStroke(ctx, PAL.green, 2, 6);

      // 3. Axión: onda dorada a(x)
      ctx.beginPath();
      for (let u = -42; u <= 42; u += 2) {
        const px = x3 + u;
        const py = y - 20 * Math.sin(u * 0.14 - ph * 1.6) * Math.exp(-u * u / 1800);
        u === -42 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
      }
      glowStroke(ctx, PAL.yellow, 2.2, 8);

      // Flechas de acoplo estilo 3B1B
      arrow(ctx, x1 + 46, y, x2 - 38, y, PAL.blue, 6, 1.6);
      arrow(ctx, x2 - 38, y, x1 + 46, y, PAL.blue, 6, 1.6);
      arrow(ctx, x2 + 38, y, x3 - 48, y, PAL.yellow, 6, 1.8);
      arrow(ctx, x3 - 48, y, x2 + 38, y, PAL.yellow, 6, 1.8);

      ctx.font = '600 12.5px "IBM Plex Sans", sans-serif';
      ctx.textAlign = 'center';
      ctx.fillStyle = PAL.blue; ctx.fillText('Gravedad g_μν', x1, 26);
      ctx.fillStyle = PAL.green; ctx.fillText('Fluido ε₀(n)', x2, 26);
      ctx.fillStyle = PAL.yellow; ctx.fillText('Axión a(x)', x3, 26);

      ctx.font = '500 11.5px "IBM Plex Mono", monospace';
      ctx.fillStyle = PAL.muted;
      ctx.fillText('T_μν', (x1 + x2) / 2, y - 12);
      ctx.fillStyle = PAL.yellow;
      ctx.fillText('n [m*(a) − m_N]', (x2 + x3) / 2, y - 12);
      ctx.fillStyle = PAL.muted;
      ctx.fillText('-R / 16πG', x1, H - 16);
      ctx.fillText('conservación bariónica', x2, H - 16);
      ctx.fillText('½(∂a)² − V(a)', x3, H - 16);
    }

    // Escena 1: Fondo (V_eff invierte su mínimo a θ=π + perfil radial θ₀(r))
    function scene1(ph) {
      const mid = W * 0.48;
      // Izquierda: Potencial efectivo V_eff(θ)
      const lx0 = 36, lx1 = mid - 24, ly0 = H - 34, ly1 = 42, lw = lx1 - lx0, lh = ly0 - ly1;
      arrow(ctx, lx0, ly0, lx1 + 8, ly0, PAL.muted, 5, 1.2);
      arrow(ctx, lx0 + lw / 2, ly0 + 6, lx0 + lw / 2, ly1 - 6, PAL.muted, 5, 1.2);

      // curva exterior (n < n_c, mínimo en θ=0) punteada azul
      ctx.save();
      ctx.setLineDash([4, 4]);
      ctx.beginPath();
      for (let i = 0; i <= 80; i++) {
        const th = (i / 80 * 2 - 1) * Math.PI;
        const v = 0.55 * (1 - Math.cos(th));
        const px = lx0 + (i / 80) * lw, py = ly0 - v * lh * 0.75 - 6;
        i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
      }
      ctx.strokeStyle = PAL.blue; ctx.lineWidth = 1.6; ctx.stroke();
      ctx.restore();

      // curva interior (n > n_c, mínimo en θ=±π) dorada
      ctx.beginPath();
      for (let i = 0; i <= 80; i++) {
        const th = (i / 80 * 2 - 1) * Math.PI;
        const v = 0.12 + 0.68 * (1 + Math.cos(th)) / 2;
        const px = lx0 + (i / 80) * lw, py = ly0 - v * lh * 0.78;
        i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
      }
      glowStroke(ctx, PAL.yellow, 2.2, 6);

      // Punto brillante en el mínimo θ = π
      const bx = lx1 - 2 + 2 * Math.sin(ph * 2), by = ly0 - 0.12 * lh * 0.78;
      ctx.fillStyle = PAL.yellow;
      ctx.beginPath(); ctx.arc(bx, by, 5, 0, 2 * Math.PI); ctx.fill();

      ctx.font = '600 12px "IBM Plex Sans", sans-serif';
      ctx.textAlign = 'left';
      ctx.fillStyle = PAL.yellow; ctx.fillText('Interior (n > n_c): mín. en θ = π', lx0, 24);
      ctx.font = '500 11px "IBM Plex Mono", monospace';
      ctx.fillStyle = PAL.blue; ctx.fillText('Fuera: mín. en θ = 0', lx0, 40);
      ctx.textAlign = 'center'; ctx.fillStyle = PAL.muted;
      ctx.fillText('−π', lx0 + 4, ly0 + 16);
      ctx.fillText('0', lx0 + lw / 2, ly0 + 16);
      ctx.fillStyle = PAL.yellow; ctx.fillText('+π', lx1 - 2, ly0 + 16);

      // Derecha: Perfil radial real θ₀(r)/π desde cw.prof
      const rx0 = mid + 24, rx1 = W - 26, rw = rx1 - rx0;
      const prof = (DATA.axcow && DATA.axcow.prof) ? DATA.axcow.prof[5] : null;
      const rSurf = rx0 + (1 / 2.2) * rw;

      // Banda de la capa de transición
      ctx.fillStyle = 'rgba(247, 217, 111, 0.14)';
      ctx.fillRect(rSurf - rw * 0.11, ly1, rw * 0.14, lh);
      ctx.save();
      ctx.setLineDash([3, 4]);
      ctx.strokeStyle = 'rgba(230,237,246,0.35)';
      ctx.beginPath(); ctx.moveTo(rSurf, ly1); ctx.lineTo(rSurf, ly0); ctx.stroke();
      ctx.restore();

      arrow(ctx, rx0, ly0, rx1 + 8, ly0, PAL.muted, 5, 1.2);
      arrow(ctx, rx0, ly0, rx0, ly1 - 6, PAL.muted, 5, 1.2);

      ctx.beginPath();
      const N = 44;
      for (let i = 0; i <= N; i++) {
        const u = i / N;
        const val = prof ? prof[i] : 1 / (1 + Math.exp((u - 0.42) * 18));
        const px = rx0 + u * rw, py = ly0 - val * (lh - 12);
        i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
      }
      glowStroke(ctx, PAL.yellow, 2.4, 7);

      ctx.font = '600 12px "IBM Plex Sans", sans-serif';
      ctx.textAlign = 'left'; ctx.fillStyle = PAL.ink;
      ctx.fillText('Perfil estático θ₀(r) (Ec. 2.6d)', rx0, 24);
      ctx.font = '500 11px "IBM Plex Mono", monospace';
      ctx.fillStyle = PAL.yellow;
      ctx.fillText('θ₀ ≈ π (condensado)', rx0 + 8, ly1 + 22);
      ctx.fillText('capa', rSurf - 16, ly1 + 12);
      ctx.fillStyle = PAL.muted;
      ctx.fillText('r = R', rSurf - 14, ly0 + 16);
      ctx.fillText('r →', rx1 - 14, ly0 + 16);
    }

    // Escena 2: Cowling (El fluido δn_b agita al axión solo en la capa donde ∂m*/∂θ ≠ 0)
    function scene2(ph) {
      const x0 = 42, x1 = W - 32, yMid = H * 0.56, w = x1 - x0;
      const xLayer = x0 + w * 0.44, wLayer = w * 0.11;

      // Banda dorada de acoplo ∂m*/∂θ ∝ sin θ₀
      const g = ctx.createLinearGradient(xLayer - wLayer, 0, xLayer + wLayer, 0);
      g.addColorStop(0, 'rgba(247, 217, 111, 0)');
      g.addColorStop(0.5, 'rgba(247, 217, 111, 0.25)');
      g.addColorStop(1, 'rgba(247, 217, 111, 0)');
      ctx.fillStyle = g;
      ctx.fillRect(xLayer - wLayer, 32, wLayer * 2, H - 56);

      // Eje radial
      arrow(ctx, x0, yMid, x1 + 8, yMid, PAL.muted, 5, 1.2);

      // 1. Onda del fluido δn_b(r) dentro de la estrella (verde Manim)
      ctx.beginPath();
      const iMax = Math.floor(w * 0.48);
      for (let u = 0; u <= iMax; u += 2) {
        const frac = u / iMax;
        const y = yMid - 34 * Math.sin(frac * 2 * Math.PI) * Math.cos(ph * 1.8);
        u === 0 ? ctx.moveTo(x0 + u, y) : ctx.lineTo(x0 + u, y);
      }
      glowStroke(ctx, PAL.green, 2.0, 5);

      // 2. Campana de acoplo ∂m*/∂θ (punteada amarilla)
      ctx.save();
      ctx.setLineDash([3, 3]);
      ctx.beginPath();
      for (let u = 0; u <= w; u += 3) {
        const px = x0 + u;
        const bell = Math.exp(-((px - xLayer) ** 2) / ((wLayer * 0.55) ** 2));
        const py = yMid - bell * 54;
        u === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
      }
      ctx.strokeStyle = PAL.yellow; ctx.lineWidth = 1.5; ctx.stroke();
      ctx.restore();

      // 3. Onda axiónica δθ(r) emitida hacia r > R (dorado/coral)
      ctx.beginPath();
      for (let u = Math.floor(w * 0.36); u <= w; u += 2) {
        const px = x0 + u;
        const env = 1 / (1 + Math.exp(-(px - (xLayer - 14)) * 0.18));
        const decay = Math.min(1, 90 / Math.max(60, px - xLayer));
        const py = yMid - 38 * env * decay * Math.cos((px - xLayer) * 0.14 - ph * 1.8);
        u === Math.floor(w * 0.36) ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
      }
      glowStroke(ctx, PAL.yellow, 2.3, 7);

      ctx.font = '600 12px "IBM Plex Sans", sans-serif';
      ctx.textAlign = 'left';
      ctx.fillStyle = PAL.green; ctx.fillText('Interior: fluido oscila (δn_b)', x0, 24);
      ctx.textAlign = 'center';
      ctx.fillStyle = PAL.yellow; ctx.fillText('Capa: ∂m*/∂θ ∝ sin θ₀ ≠ 0', xLayer, 24);
      ctx.textAlign = 'right';
      ctx.fillStyle = PAL.yellow; ctx.fillText('Exterior: onda δθ escapa si ω > m_a →', x1, 24);

      ctx.font = '500 11px "IBM Plex Mono", monospace';
      ctx.textAlign = 'center'; ctx.fillStyle = PAL.muted;
      ctx.fillText('centro (r=0)', x0 + 28, H - 12);
      ctx.fillText('superficie (r=R)', xLayer + 12, H - 12);
      ctx.fillText('espacio exterior (r > R)', x0 + w * 0.78, H - 12);
    }

    // Escena 3: Más allá de Cowling (Comparación directa: Cowling vs Métrica dinámica)
    function scene3(ph) {
      const cx1 = W * 0.26, cx2 = W * 0.74, cy = H * 0.56, r0 = 30;
      const cosw = Math.cos(ph * 1.8);

      for (const [cx, dyn, title, sub, col] of [
        [cx1, false, 'Cowling (δg_μν = 0)', 'estrella vibra · métrica quieta', PAL.muted],
        [cx2, true, 'Métrica dinámica (H₀, H₁, K)', 'el espacio-tiempo vibra y emite GW', PAL.blue]
      ]) {
        // Malla polar/cartesiana alrededor de la estrella
        ctx.strokeStyle = dyn ? 'rgba(88, 196, 221, 0.28)' : 'rgba(139, 151, 171, 0.16)';
        ctx.lineWidth = 1;
        for (let gx = -3; gx <= 3; gx++) {
          ctx.beginPath();
          for (let gy = -20; gy <= 20; gy += 2) {
            const dx = gx * 22, dy = gy * 3.4, rho = Math.hypot(dx, dy);
            let wx = cx + dx, wy = cy + dy;
            if (dyn && rho > r0 * 0.9) {
              const c2 = (dy * dy) / (rho * rho + 1);
              const h = 5.2 * (1 - c2) * Math.cos(ph * 1.8 - (rho - r0) * 0.16);
              wx += h * (dx / (rho + 1)); wy -= h * (dy / (rho + 1));
            }
            gy === -20 ? ctx.moveTo(wx, wy) : ctx.lineTo(wx, wy);
          }
          ctx.stroke();
        }
        // Frentes de onda gravitacional en el lado dinámico
        if (dyn) {
          for (let ring = 0; ring < 3; ring++) {
            const rr = r0 + 10 + ((ph * 18 + ring * 26) % 78);
            const alpha = Math.max(0, 1 - (rr - r0) / 82) * 0.65;
            ctx.strokeStyle = `rgba(88, 196, 221, ${alpha.toFixed(2)})`;
            ctx.lineWidth = 1.8;
            ctx.beginPath();
            ctx.ellipse(cx, cy, rr * 1.12, rr * 0.88, 0, 0, 2 * Math.PI);
            ctx.stroke();
          }
        }
        // Estrella cuadrupolar
        ctx.beginPath();
        for (let a = 0; a <= 60; a++) {
          const ang = (a / 60) * 2 * Math.PI, c = Math.cos(ang), P2 = 1.5 * c * c - 0.5;
          const rr = r0 * (1 + 0.14 * P2 * cosw);
          const px = cx + rr * Math.sin(ang), py = cy - rr * c;
          a === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
        }
        ctx.fillStyle = 'rgba(28, 43, 74, 0.9)'; ctx.fill();
        glowStroke(ctx, dyn ? PAL.blue : PAL.ink, 2, 5);

        ctx.textAlign = 'center';
        ctx.font = '600 12.5px "IBM Plex Sans", sans-serif';
        ctx.fillStyle = col; ctx.fillText(title, cx, 24);
        ctx.font = '400 11.5px "IBM Plex Sans", sans-serif';
        ctx.fillStyle = PAL.muted; ctx.fillText(sub, cx, H - 12);
      }
    }

    // Escena 4: Sistema extendido (Diagrama de acoplo Lindblom-Detweiler + términos nuevos en amarillo)
    function scene4(ph) {
      const yBox = H * 0.48, bw = Math.min(155, W * 0.25), bh = 56;
      const x1 = W * 0.19, x2 = W * 0.5, x3 = W * 0.81;

      function box(cx, title, vars, col) {
        ctx.save();
        ctx.fillStyle = 'rgba(20, 27, 39, 0.95)';
        ctx.strokeStyle = col;
        ctx.lineWidth = 1.8;
        ctx.shadowColor = col; ctx.shadowBlur = 6;
        ctx.beginPath();
        ctx.roundRect ? ctx.roundRect(cx - bw / 2, yBox - bh / 2, bw, bh, 7) : ctx.rect(cx - bw / 2, yBox - bh / 2, bw, bh);
        ctx.fill(); ctx.stroke();
        ctx.restore();
        ctx.textAlign = 'center';
        ctx.font = '600 12.5px "IBM Plex Sans", sans-serif';
        ctx.fillStyle = col; ctx.fillText(title, cx, yBox - 6);
        ctx.font = '500 11px "IBM Plex Mono", monospace';
        ctx.fillStyle = PAL.ink; ctx.fillText(vars, cx, yBox + 14);
      }

      // Flecha clásica LD (azul) entre Métrica y Fluido
      arrow(ctx, x1 + bw / 2 + 4, yBox - 8, x2 - bw / 2 - 4, yBox - 8, PAL.blue, 5.5, 1.6);
      arrow(ctx, x2 - bw / 2 - 4, yBox + 8, x1 + bw / 2 + 4, yBox + 8, PAL.blue, 5.5, 1.6);

      // Flechas de términos nuevos (amarillo brillante) con pulso viajero
      arrow(ctx, x2 + bw / 2 + 4, yBox - 8, x3 - bw / 2 - 4, yBox - 8, PAL.yellow, 5.5, 2.0);
      arrow(ctx, x3 - bw / 2 - 4, yBox + 8, x2 + bw / 2 + 4, yBox + 8, PAL.yellow, 5.5, 2.0);

      // Arco directo Métrica <-> Axión por debajo (flujo 8πG r f_a a_0' δθ y término H_0)
      ctx.save();
      ctx.strokeStyle = PAL.yellow;
      ctx.lineWidth = 2;
      ctx.shadowColor = PAL.yellow; ctx.shadowBlur = 7;
      ctx.beginPath();
      ctx.moveTo(x1, yBox + bh / 2 + 2);
      ctx.quadraticCurveTo(x2, H - 8, x3, yBox + bh / 2 + 2);
      ctx.stroke();
      ctx.restore();

      // Pulso luminoso recorriendo el arco nuevo Métrica <-> Axión
      const u = 0.5 + 0.45 * Math.sin(ph * 1.6);
      const px = (1 - u) ** 2 * x1 + 2 * (1 - u) * u * x2 + u ** 2 * x3;
      const py = (1 - u) ** 2 * (yBox + bh / 2 + 2) + 2 * (1 - u) * u * (H - 8) + u ** 2 * (yBox + bh / 2 + 2);
      ctx.fillStyle = PAL.yellow;
      ctx.beginPath(); ctx.arc(px, py, 4.5, 0, 2 * Math.PI); ctx.fill();

      box(x1, 'Métrica (LD)', 'H₁, K, H₀', PAL.blue);
      box(x2, 'Fluido (LD)', 'W, X, V', PAL.green);
      box(x3, 'Axión (KG)', 'δθ, δθ′', PAL.yellow);

      ctx.font = '500 11px "IBM Plex Mono", monospace';
      ctx.textAlign = 'center';
      ctx.fillStyle = PAL.blue;
      ctx.fillText('Lindblom–Detweiler (1983)', (x1 + x2) / 2, 26);
      ctx.fillStyle = PAL.yellow;
      ctx.fillText('n_b M′ Δθ (nuevo)', (x2 + x3) / 2, 26);
      ctx.fillText('Acoplo directo métrica–axión: +8πG r fₐ a₀′ δθ  y  −2e⁻²λ Φ′ θ₀′ H₀', x2, H - 6);
    }

    // Escena 5: Modos cuasinormales (Plano complejo ω + dos canales de amortiguamiento)
    function scene5(ph) {
      const mid = W * 0.48;
      // Izquierda: Plano complejo (Re ν, Im ν) con umbral ν_a
      const lx0 = 40, lx1 = mid - 20, lyTop = 36, lyBot = H - 28;
      arrow(ctx, lx0, lyTop + 14, lx1 + 6, lyTop + 14, PAL.muted, 5, 1.3);
      arrow(ctx, lx0 + 12, lyTop, lx0 + 12, lyBot + 6, PAL.muted, 5, 1.3);

      const xThresh = lx0 + (lx1 - lx0) * 0.48;
      ctx.fillStyle = 'rgba(247, 217, 111, 0.09)';
      ctx.fillRect(xThresh, lyTop + 14, lx1 - xThresh, lyBot - lyTop - 14);
      ctx.save();
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = PAL.yellow;
      ctx.beginPath(); ctx.moveTo(xThresh, lyTop + 4); ctx.lineTo(xThresh, lyBot); ctx.stroke();
      ctx.restore();

      const pts = [
        { x: lx0 + (lx1 - lx0) * 0.26, y: lyTop + 42, col: PAL.blue, lab: 'f (solo GW)' },
        { x: lx0 + (lx1 - lx0) * 0.66, y: lyTop + 68, col: PAL.yellow, lab: 'p₁ (GW + axión)' },
        { x: lx0 + (lx1 - lx0) * 0.86, y: lyTop + 108, col: PAL.red, lab: 'p₂ (τ ≈ 12 ms)' },
      ];
      for (const p of pts) {
        ctx.fillStyle = p.col;
        ctx.beginPath(); ctx.arc(p.x, p.y, 4.8 + 1.2 * Math.sin(ph * 2), 0, 2 * Math.PI); ctx.fill();
        ctx.font = '500 11px "IBM Plex Mono", monospace';
        ctx.textAlign = 'left';
        ctx.fillText(p.lab, p.x - 18, p.y - 9);
      }
      ctx.font = '600 12px "IBM Plex Sans", sans-serif';
      ctx.fillStyle = PAL.ink; ctx.fillText('Plano complejo: ω = 2π(Re ν − i |Im ν|)', lx0, 22);
      ctx.font = '500 11px "IBM Plex Mono", monospace';
      ctx.fillStyle = PAL.yellow; ctx.fillText('umbral Re ν = νₐ', xThresh + 4, lyBot - 4);
      ctx.fillStyle = PAL.muted; ctx.fillText('−Im ν ↓', lx0 + 18, lyBot - 4);

      // Derecha: Onda amortiguada e^(-t/τ) cos(ωt)
      const rx0 = mid + 24, rx1 = W - 24, ryMid = H * 0.54, rw = rx1 - rx0;
      arrow(ctx, rx0, ryMid, rx1 + 6, ryMid, PAL.muted, 5, 1.2);
      ctx.save();
      ctx.setLineDash([3, 3]);
      ctx.strokeStyle = 'rgba(252, 98, 85, 0.6)';
      for (const sgn of [-1, 1]) {
        ctx.beginPath();
        for (let u = 0; u <= rw; u += 3) {
          const env = 54 * Math.exp(-2.1 * (u / rw));
          u === 0 ? ctx.moveTo(rx0 + u, ryMid + sgn * env) : ctx.lineTo(rx0 + u, ryMid + sgn * env);
        }
        ctx.stroke();
      }
      ctx.restore();

      ctx.beginPath();
      for (let u = 0; u <= rw; u += 1.5) {
        const frac = u / rw;
        const y = ryMid - 54 * Math.exp(-2.1 * frac) * Math.cos(frac * 7.5 * Math.PI - ph * 2.2);
        u === 0 ? ctx.moveTo(rx0 + u, y) : ctx.lineTo(rx0 + u, y);
      }
      glowStroke(ctx, PAL.blue, 2.2, 6);

      ctx.font = '600 12px "IBM Plex Sans", sans-serif';
      ctx.textAlign = 'left'; ctx.fillStyle = PAL.ink;
      ctx.fillText('Modo cuasinormal: vibración que se apaga', rx0, 22);
      ctx.font = '500 11px "IBM Plex Mono", monospace';
      ctx.fillStyle = PAL.red;
      ctx.fillText('envolvente ±e^(−t/τ),  1/τ = 2π |Im ν|', rx0 + 6, H - 14);
    }

    const scenes = [scene0, scene1, scene2, scene3, scene4, scene5];

    function draw() {
      if (!ctx) return;
      drawBackgroundGrid();
      ctx.save();
      ctx.globalAlpha = Math.min(1, stepAnim);
      scenes[k](t);
      ctx.restore();
    }

    function frame(ts) {
      if (last !== null) {
        const dt = Math.min(0.05, (ts - last) / 1000);
        if (playing) t += dt;
        if (stepAnim < 1) stepAnim = Math.min(1, stepAnim + dt * 3.2);
      }
      last = ts;
      if (visible) draw();
      if ((playing || stepAnim < 1) && visible) requestAnimationFrame(frame); else last = null;
    }

    function show(n) {
      k = (n + steps.length) % steps.length;
      stepAnim = 0.15;
      steps.forEach((s, i) => { s.classList.toggle('on', i === k); s.setAttribute('aria-hidden', String(i !== k)); });
      dots.forEach((d, i) => { d.setAttribute('aria-current', i === k ? 'step' : 'false'); d.classList.toggle('done', i < k); });
      draw();
      if (last === null && visible) requestAnimationFrame(frame);
    }

    prev.addEventListener('click', () => show(k - 1));
    next.addEventListener('click', () => show(k + 1));
    dots.forEach((d, i) => d.addEventListener('click', () => show(i)));
    if (playBtn) {
      playBtn.addEventListener('click', () => {
        playing = !playing;
        playBtn.textContent = playing ? 'Pausa' : 'Animar';
        if (playing && last === null) requestAnimationFrame(frame);
      });
    }
    root.addEventListener('keydown', e => { if (e.key === 'ArrowRight') show(k + 1); if (e.key === 'ArrowLeft') show(k - 1); });
    if (cv) {
      new IntersectionObserver(es => { visible = es[0].isIntersecting; if (visible && (playing || stepAnim < 1) && last === null) requestAnimationFrame(frame); }).observe(cv);
      let lastW = 0;
      new ResizeObserver(() => {
        const w = cv.parentElement.clientWidth;
        if (Math.abs(w - lastW) > 1) { lastW = w; layout(); draw(); }
      }).observe(cv.parentElement);
      layout();
    }
    show(0);
  }

  document.querySelectorAll('[data-anim=mode]').forEach(modeAnim);
  document.querySelectorAll('[data-anim=chain]').forEach(stepper);
})();
