// Animaciones de la página: la anatomía de un modo (datos de los pasos 9 y 10 con estética 3Blue1Brown / Manim)
// y la cadena de seis ecuaciones en LaTeX. Sin bibliotecas externas: Canvas 2D y CSS.
(() => {
  const DATA = JSON.parse(document.getElementById('figdata').textContent);
  const REDUCE = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const PAL = {
    bg: [14, 18, 28], star: [28, 43, 74], cond: [247, 217, 111],
    gwp: [88, 196, 221], gwn: [154, 114, 172],
    axp: [247, 217, 111], axn: [252, 98, 85],
    blue: '#58c4dd', yellow: '#f7d96f', red: '#fc6255',
    ink: '#e6edf6', muted: '#8b97ab'
  };

  function arrow(ctx, x1, y1, x2, y2, color, size = 5.5, lw = 1.4) {
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

  function glowStroke(ctx, color, width, glow = 6) {
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.shadowColor = color;
    ctx.shadowBlur = glow;
    ctx.stroke();
    ctx.shadowBlur = 0;
  }

  // ------------------------------------------------------------------ anatomía de un modo
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
      H = Math.round(Math.min(445, Math.max(305, W * 0.52)));
      cv.width = W * dpr; cv.height = H * dpr; cv.style.height = H + 'px';
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      R = Math.min(H * 0.235, W * 0.13); cx = W / 2; cy = H / 2 - 4;
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

    // perfil radial del desplazamiento (estilizado): f sin nodos, p1 uno, p2 dos
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
      const amp = 0.072, cosw = Math.cos(ph), sinw = Math.sin(ph);
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
            d[o + q] = mix(mix(PAL.star[q], PAL.cond[q], th ** 1.5 * 0.75), 255, 0.1 * disp) * (0.7 + 0.3 * tone);
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
          col = v > 0 ? PAL.axp : PAL.axn; a = Math.min(1, Math.abs(v)) * 0.8;
        }
        for (let q = 0; q < 3; q++) d[o + q] = col ? mix(c[q], col[q], a) : c[q];
        d[o + 3] = 255;
      }
      octx.putImageData(img, 0, 0);
      ctx.imageSmoothingEnabled = true;
      ctx.drawImage(off, 0, 0, off.width * S, off.height * S);

      // Retícula de espacio-tiempo 3Blue1Brown en el hemisferio izquierdo (oscila con h+ en métrica dinámica, rígida en Cowling)
      ctx.save();
      ctx.strokeStyle = st.dyn ? 'rgba(88, 196, 221, 0.17)' : 'rgba(139, 151, 171, 0.10)';
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
            const h = 4.2 * (1 - c2) * Math.cos(ph - 2 * Math.PI * (rho - R) / lamGW) * Math.min(1, 1.5 * R / rho);
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
            const h = 4.2 * (1 - c2) * Math.cos(ph - 2 * Math.PI * (rho - R) / lamGW) * Math.min(1, 1.5 * R / rho);
            wx += h * (dx / rho);
            wy -= h * (dy / rho);
          }
          gx === -32 ? ctx.moveTo(wx, wy) : ctx.lineTo(wx, wy);
        }
        ctx.stroke();
      }
      ctx.restore();

      // Capas del fluido + superficie
      for (const r0 of [0.25, 0.5, 0.75, 1]) {
        ctx.beginPath();
        for (let a2 = 0; a2 <= 72; a2++) {
          const t = a2 / 72 * 2 * Math.PI, c = Math.cos(t), P2 = 1.5 * c * c - 0.5;
          const rr = R * (r0 + amp * xi[m.i](r0) * P2 * cosw);
          const X = cx + rr * Math.sin(t), Y = cy - rr * c;
          a2 ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y);
        }
        if (r0 === 1) {
          glowStroke(ctx, 'rgba(230,237,246,0.95)', 2.1, 5);
        } else {
          ctx.strokeStyle = 'rgba(14,18,28,0.45)';
          ctx.lineWidth = 1;
          ctx.stroke();
        }
      }

      // Nodos radiales de p1 y p2 en azul Manim (#58c4dd)
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

      // Vectores de velocidad cuadrupolar del fluido dentro de la estrella
      const vphase = -sinw;
      if (Math.abs(vphase) > 0.12) {
        const rSample = m.i === 0 ? 0.62 : m.i === 1 ? 0.84 : 0.68;
        const sgn = Math.sign(xi[m.i](rSample)) * vphase;
        for (const ang of [0, Math.PI / 2, Math.PI, 1.5 * Math.PI]) {
          const c = Math.cos(ang), s = Math.sin(ang), P2 = 1.5 * c * c - 0.5;
          const rBase = R * rSample;
          const x1 = cx + rBase * s, y1 = cy - rBase * c;
          const dr = 10 * P2 * sgn;
          arrow(ctx, x1, y1, x1 + dr * s, y1 - dr * c, 'rgba(230,237,246,0.82)', 4.5, 1.3);
        }
      }

      // Indicador vectorial de la capa de transición (θ₀ ≈ π/2)
      const xTrans = transitionRadius(prof);
      if (xTrans !== null && W > 560) {
        const rT = R * xTrans;
        ctx.save();
        ctx.setLineDash([2, 3]);
        ctx.strokeStyle = PAL.yellow;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.arc(cx, cy, rT, -0.65, 0.65);
        ctx.stroke();
        ctx.restore();

        const tx = cx + rT * Math.cos(-0.4), ty = cy + rT * Math.sin(-0.4);
        const lx = cx + R * 1.45, ly = cy - R * 0.9;
        arrow(ctx, lx - 4, ly + 4, tx + 3, ty - 2, PAL.yellow, 5, 1.3);
        ctx.font = '500 11.5px "IBM Plex Mono", monospace';
        ctx.fillStyle = PAL.yellow;
        ctx.textAlign = 'left';
        ctx.fillText('capa de transición', lx, ly);
      }

      // Rótulos superiores
      ctx.font = '500 13px "IBM Plex Sans", system-ui, sans-serif';
      ctx.fillStyle = PAL.ink; ctx.textAlign = 'left';
      ctx.fillText('ondas gravitacionales', 14, 22);
      ctx.textAlign = 'right';
      ctx.fillText('radiación axiónica, δθ', W - 14, 22);
      ctx.font = '400 12px "IBM Plex Sans", system-ui, sans-serif';
      ctx.fillStyle = PAL.muted;
      ctx.textAlign = 'left';
      ctx.fillText(st.dyn ? 'h₊ ∝ sin²ϑ' : 'Cowling: métrica fija, no hay', 14, 40);
      ctx.textAlign = 'right';
      ctx.fillText(m.nua === null ? 'sin axión' : m.radiates ? 'ν > νₐ: se propaga' : 'ν < νₐ: ligada, decae exponencialmente', W - 14, 40);

      // Traza de ringdown e^(-t/τ) cos(ωt) en esquina inferior izquierda
      if (W > 580) {
        const bx = 14, by = H - 64, bw = 168, bh = 40;
        ctx.save();
        ctx.fillStyle = 'rgba(20, 27, 39, 0.85)';
        ctx.strokeStyle = '#2c3850';
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.roundRect ? ctx.roundRect(bx, by, bw, bh, 5) : ctx.rect(bx, by, bw, bh);
        ctx.fill(); ctx.stroke();

        const decay = m.im > 1e-4 ? Math.min(2.6, 0.32 * Math.log10(1 + m.im * 18)) : 0;
        const midY = by + bh / 2, ampY = bh * 0.35;

        if (decay > 0.01) {
          ctx.setLineDash([2, 2]);
          ctx.strokeStyle = 'rgba(252, 98, 85, 0.5)';
          for (const sgn of [-1, 1]) {
            ctx.beginPath();
            for (let u = 0; u <= bw - 14; u += 3) {
              const env = Math.exp(-decay * (u / (bw - 14)));
              const px = bx + 7 + u, py = midY + sgn * ampY * env;
              u === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py);
            }
            ctx.stroke();
          }
          ctx.setLineDash([]);
        }

        ctx.beginPath();
        for (let u = 0; u <= bw - 14; u += 1.5) {
          const frac = u / (bw - 14);
          const env = Math.exp(-decay * frac);
          const y = midY - ampY * env * Math.cos(frac * 6.5 * Math.PI - ph);
          u === 0 ? ctx.moveTo(bx + 7 + u, y) : ctx.lineTo(bx + 7 + u, y);
        }
        glowStroke(ctx, m.radiates ? PAL.yellow : PAL.blue, 1.5, 4);
        ctx.restore();
      }

      ctx.textAlign = 'center';
      ctx.font = '400 12px "IBM Plex Sans", system-ui, sans-serif';
      ctx.fillStyle = PAL.muted;
      if (W > 740) ctx.fillText('fuera de la estrella las distancias están comprimidas; la proporción entre longitudes de onda es la real', W / 2 + 50, H - 12);
    }

    function readout() {
      const m = current(), lab = ['f', 'p₁', 'p₂'][m.i];
      ui.nu.textContent = m.nua === null ? 'sin axión' : `νₐ = ${m.nua} kHz`;
      ui.freq.textContent = `${lab}: ${m.re.toFixed(3)} kHz`;
      if (m.im > 1e-3) {
        const tau = 1 / (2 * Math.PI * m.im);
        ui.damp.textContent = `|Im ν| = ${m.im < 0.1 ? m.im.toPrecision(2) : m.im.toFixed(2)} Hz · τ = ${tau >= 1 ? tau.toFixed(1) + ' s' : (tau * 1e3).toPrecision(3) + ' ms'}`;
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

  // ------------------------------------------------------------------ cadena de ecuaciones
  function stepper(root) {
    const steps = [...root.querySelectorAll('.step')], dots = [...root.querySelectorAll('.rail button')];
    const prev = root.querySelector('[data-prev]'), next = root.querySelector('[data-next]');
    let k = 0;
    function show(n) {
      k = (n + steps.length) % steps.length;
      steps.forEach((s, i) => { s.classList.toggle('on', i === k); s.setAttribute('aria-hidden', String(i !== k)); });
      dots.forEach((d, i) => { d.setAttribute('aria-current', i === k ? 'step' : 'false'); d.classList.toggle('done', i < k); });
    }
    prev.addEventListener('click', () => show(k - 1));
    next.addEventListener('click', () => show(k + 1));
    dots.forEach((d, i) => d.addEventListener('click', () => show(i)));
    root.addEventListener('keydown', e => { if (e.key === 'ArrowRight') show(k + 1); if (e.key === 'ArrowLeft') show(k - 1); });
    show(0);
  }

  document.querySelectorAll('[data-anim=mode]').forEach(modeAnim);
  document.querySelectorAll('[data-anim=chain]').forEach(stepper);
})();
