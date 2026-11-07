"""Genera report/informe_magnetar.html y results/informe_magnetar.html desde checks_magnetar.json y fig_magnetar.json."""
import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / 'results'
REP = ROOT / 'report'


def build():
    checks_data = json.loads((RES / 'checks_magnetar.json').read_text())
    fig_data = json.loads((RES / 'fig_magnetar.json').read_text())

    checks_html = []
    for c in checks_data['checks']:
        kind = c.get('kind', 'symbolic')
        if not c['passed']:
            chip = '<span class="chip chip-bad">no cuadra</span>'
        elif kind == 'negative_control':
            chip = '<span class="chip chip-obs">control neg.</span>'
        else:
            chip = '<span class="chip chip-ok">coincide</span>'
        checks_html.append(
            f'<li>{chip}<span class="ck">{html.escape(c["name"])}</span></li>'
        )
    checks_list_str = '\n'.join(checks_html)

    rows_html = []
    for r in fig_data['rows']:
        nu = r['nu_a_khz']
        b14 = r['magnetic_fields']['1e+14']
        b514 = r['magnetic_fields']['5e+14']
        b15 = r['magnetic_fields']['1e+15']
        rows_html.append(
            f'<tr>'
            f'<td class="num"><b>{nu:.1f}</b></td>'
            f'<td class="num">{r["M_msun"]:.3f}</td>'
            f'<td class="num">{r["R_km"]:.2f}</td>'
            f'<td class="num">{r["r_wall_km"]:.2f}</td>'
            f'<td class="num">{r["ws_amplification"]:.4f} <span class="muted">({r["delta_ws_vs_no_axion_pct"]:+.2f}%)</span></td>'
            f'<td class="num">{r["k_crit_max_cm_inv"]:.2f}</td>'
            f'<td class="num">{r["lambda_crit_cm"]:.2f}</td>'
            f'<td class="num">{r["mu5_eff_max_ev"]:.4f}</td>'
            f'<td class="num">{b15["delta_nu_hz"]:+.1f}</td>'
            f'<td class="num">{b14["L_gamma_erg_s"]:.2e}</td>'
            f'<td class="num">{b514["L_gamma_erg_s"]:.2e}</td>'
            f'<td class="num"><b>{b15["L_gamma_erg_s"]:.2e}</b></td>'
            f'</tr>'
        )
    table_rows_str = '\n'.join(rows_html)
    fig_json_inline = json.dumps(fig_data, separators=(',', ':')).replace('</', '<\\/')

    page = f"""<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Unificación Axión + Magnetar: Chern–Simons, corriente quiral y conversión resonante</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {{
  --paper: #f4f6f9;
  --sheet: #ffffff;
  --ink: #161e2e;
  --muted: #586378;
  --rule: #d8dde6;
  --accent: #2b4bb3;
  --ok: #1d7451;
  --ok-bg: #e4f3eb;
  --flag: #9a5400;
  --flag-bg: #fdf0dc;
  --bad: #a3262f;
  --bad-bg: #fbe5e6;
  --code-bg: #eaeef4;
  --c1: #2b4bb3; --c2: #c2410c; --c3: #7c3aed; --c4: #0d9488;
  --f-display: "Newsreader", "Iowan Old Style", Georgia, serif;
  --f-body: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  --f-mono: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper: #0d121b; --sheet: #141b27; --ink: #e3e8f1; --muted: #98a3b7; --rule: #263042;
    --accent: #93a8ff; --ok: #62c99c; --ok-bg: #12291f; --flag: #f1a84d; --flag-bg: #2b2112;
    --bad: #f08a90; --bad-bg: #2d1518; --code-bg: #1b2433; color-scheme: dark;
    --c1: #6283ee; --c2: #d9733a; --c3: #a174ec; --c4: #2dd4bf;
  }}
}}
:root[data-theme="dark"] {{
  --paper: #0d121b; --sheet: #141b27; --ink: #e3e8f1; --muted: #98a3b7; --rule: #263042;
  --accent: #93a8ff; --ok: #62c99c; --ok-bg: #12291f; --flag: #f1a84d; --flag-bg: #2b2112;
  --bad: #f08a90; --bad-bg: #2d1518; --code-bg: #1b2433; color-scheme: dark;
  --c1: #6283ee; --c2: #d9733a; --c3: #a174ec; --c4: #2dd4bf;
}}
* {{ box-sizing: border-box; }}
html {{ scroll-behavior: smooth; }}
body {{
  margin: 0;
  background: var(--paper); color: var(--ink);
  font: 400 1rem/1.62 var(--f-body);
  padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px));
  padding-block: 2rem 4rem;
}}
.wrap {{ display: grid; gap: 2.8rem 0;
  grid-template-columns: [full-start] minmax(0, 1fr) [content-start] minmax(0, 70ch) [content-end] minmax(0, 1fr) [full-end]; }}
.wrap > * {{ grid-column: content; min-width: 0; }}
.wrap > .wide {{ grid-column: full; width: min(1080px, 100%); justify-self: center; }}
h1, h2, h3 {{ font-family: var(--f-display); font-weight: 600; text-wrap: balance; line-height: 1.2; margin: 0; }}
h1 {{ font-size: clamp(1.95rem, 4.6vw, 2.8rem); letter-spacing: -0.01em; }}
h2 {{ font-size: 1.58rem; }}
h3 {{ font-size: 1.15rem; }}
p {{ margin: 0; }}
a {{ color: var(--accent); text-underline-offset: 2px; }}
.eyebrow {{ font: 500 0.78rem/1.4 var(--f-mono); letter-spacing: 0.06em; text-transform: uppercase; color: var(--muted); }}
.stack {{ display: grid; gap: 1rem; }}
.lede {{ font-size: 1.12rem; color: var(--ink); }}
.muted {{ color: var(--muted); }}
section[id], h3[id], footer[id] {{ scroll-margin-top: 4.2rem; }}
mjx-container[display="true"] {{ overflow-x: auto; overflow-y: hidden; max-width: 100%; padding-block: 0.2rem; }}

.topbar {{
  position: sticky; top: 0.6rem; z-index: 10;
  display: flex; flex-wrap: wrap; align-items: center; justify-content: space-between;
  gap: 0.5rem 1rem; padding: 0.55rem 0.9rem;
  background: var(--sheet); border: 1px solid var(--rule); border-radius: 8px;
  box-shadow: 0 4px 14px rgb(0 0 0 / 0.06); font-size: 0.86rem;
}}
.navlinks {{ display: flex; flex-wrap: wrap; gap: 0.3rem 0.8rem; align-items: center; }}
.navlinks a {{ text-decoration: none; font-weight: 500; color: var(--ink); padding: 0.22rem 0.48rem; border-radius: 5px; }}
.navlinks a:hover {{ color: var(--accent); text-decoration: underline; }}
.navlinks a[aria-current="page"] {{ background: var(--code-bg); color: var(--accent); font-weight: 600; }}
.actions {{ display: flex; flex-wrap: wrap; gap: 0.45rem; align-items: center; }}
.btn-link, .btn-theme {{
  font: 500 0.78rem/1.2 var(--f-mono); padding: 0.38rem 0.65rem;
  border-radius: 5px; border: 1px solid var(--rule); background: var(--code-bg);
  color: var(--ink); text-decoration: none; cursor: pointer;
}}
.btn-link {{ border-color: var(--accent); color: var(--accent); background: transparent; }}

.tally {{ display: flex; flex-wrap: wrap; gap: 0.6rem 2rem; padding-block: 1rem; border-block: 1px solid var(--rule); font-variant-numeric: tabular-nums; }}
.tally div {{ display: grid; gap: 0.1rem; min-width: 0; }}
.tally b {{ font: 600 1.5rem/1.1 var(--f-display); }}
.tally span {{ font-size: 0.85rem; color: var(--muted); }}

.stage {{ background: #0e121c; color: #e6edf6; border: 1px solid #232c3d; border-radius: 8px; padding: 0.95rem; display: grid; gap: 0.75rem; min-width: 0; }}
.stage-bar, .stage-foot {{ display: flex; flex-wrap: wrap; gap: 0.5rem 1rem; align-items: center; }}
.stage .btn {{ background: none; color: #aab6c9; border: 1px solid #2c3850; border-radius: 6px; padding: 0.4rem 0.75rem; font: 500 0.84rem/1.2 var(--f-body); cursor: pointer; }}
.stage .btn:hover {{ color: #fff; }}
.stage-canvas canvas {{ display: block; width: 100%; border-radius: 4px; background: #090d14; }}
.slider {{ display: grid; gap: 0.25rem; flex: 1 1 14rem; font-size: 0.85rem; color: #aab6c9; }}
.slider input {{ width: 100%; accent-color: #58c4dd; }}
.readout {{ display: flex; flex-wrap: wrap; align-items: baseline; gap: 0.25rem 1.2rem; font-size: 0.88rem; color: #aab6c9; font-variant-numeric: tabular-nums; }}
.readout b {{ font: 600 1.12rem/1.2 var(--f-body); color: #fff; }}

.cites {{ list-style: none; margin: 0; padding: 0; display: grid; gap: 0.9rem; }}
.cites li {{ display: grid; gap: 0.25rem; padding: 0.9rem 1.1rem; background: var(--sheet); border: 1px solid var(--rule); border-left: 3px solid var(--accent); border-radius: 0 6px 6px 0; }}
.cites .tag {{ font: 500 0.72rem/1 var(--f-mono); letter-spacing: 0.05em; text-transform: uppercase; color: var(--accent); }}

.ob {{ display: grid; gap: 0.7rem; padding: 1.2rem 1.3rem; background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; }}
.ob header {{ display: flex; flex-wrap: wrap; align-items: baseline; gap: 0.4rem 0.8rem; }}
.weight {{ font: 500 0.72rem/1 var(--f-mono); letter-spacing: 0.05em; text-transform: uppercase; padding: 0.3rem 0.5rem; border-radius: 3px; }}
.w-main {{ color: var(--flag); background: var(--flag-bg); }}
.w-ok {{ color: var(--ok); background: var(--ok-bg); }}
.newterm {{ color: var(--flag); font-weight: 600; }}

.charts {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 26rem), 1fr)); gap: 1rem; }}
.chart {{ position: relative; display: grid; gap: 0.55rem; align-content: start; padding: 0.95rem 0.95rem 0.6rem; background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; }}
.chart h4 {{ margin: 0; font: 600 0.98rem/1.3 var(--f-body); }}
.chart .sub {{ font-size: 0.82rem; color: var(--muted); }}
.plot {{ display: block; width: 100%; height: auto; font: 400 11px var(--f-body); }}

.checks {{ list-style: none; margin: 0; padding: 0; display: grid; }}
.checks li {{ display: grid; grid-template-columns: 7rem 1fr; gap: 0.7rem; align-items: baseline; padding-block: 0.45rem; border-top: 1px solid var(--rule); font-size: 0.93rem; }}
.chip {{ font: 500 0.7rem/1 var(--f-mono); letter-spacing: 0.04em; text-transform: uppercase; padding: 0.28rem 0.45rem; border-radius: 3px; text-align: center; }}
.chip-ok {{ color: var(--ok); background: var(--ok-bg); }}
.chip-obs {{ color: var(--flag); background: var(--flag-bg); }}
.chip-bad {{ color: var(--bad); background: var(--bad-bg); }}

.tablewrap {{ overflow-x: auto; background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; padding: 0.6rem 0.9rem; }}
.steps {{ width: 100%; border-collapse: collapse; font-size: 0.88rem; }}
.steps th, .steps td {{ text-align: left; padding: 0.45rem 0.65rem 0.45rem 0; border-bottom: 1px solid var(--rule); vertical-align: baseline; }}
.steps th {{ font-weight: 600; color: var(--muted); font-size: 0.8rem; }}
.steps .num {{ text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; padding-left: 0.8rem; }}
footer {{ font-size: 0.85rem; color: var(--muted); border-top: 1px solid var(--rule); padding-top: 1rem; display: flex; flex-wrap: wrap; justify-content: space-between; gap: 0.5rem; }}
</style>
<script>
window.MathJax = {{ tex: {{ inlineMath: [['\\\\(', '\\\\)']], displayMath: [['\\\\[', '\\\\]']] }}, svg: {{ fontCache: 'global' }} }};
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.min.js" async></script>
</head>
<body>

<main class="wrap">
  <nav class="topbar wide" aria-label="Navegación entre informes">
    <div class="navlinks">
      <a href="informe_axiones.html">Informe base (Cowling + GR)</a>
      <a href="informe_gmodes.html">Línea 1: Modos g (LIGO/ET)</a>
      <a href="informe_rotacion.html">Línea 2: Rotación y modos r</a>
      <a href="informe_magnetar.html" aria-current="page">Línea 3: Magnetar + Axión</a>
    </div>
    <div class="actions">
      <a class="btn-link" href="guia_departamento.html">Guía departamento →</a>
      <button type="button" class="btn-theme" id="themetoggle">Tema</button>
    </div>
  </nav>

  <header class="stack">
    <p class="eyebrow">Línea 3 · Unificación de las dos líneas de J. A. Pons (2024–2026) · Paso 13</p>
    <h1>Acoplo de Chern–Simons axión-fotón, corriente quiral en la pared de dominio y oscilaciones de magnetares</h1>
    <p class="lede">Unimos los artículos de condensados de axiones ligeros de Gómez-Bañón, Pnigouras y Pons (<a href="https://arxiv.org/abs/2408.13247">PRL 2024</a>, <a href="https://arxiv.org/abs/2510.16245">arXiv:2510.16245</a>, <a href="https://arxiv.org/abs/2610.37308">arXiv:2610.37308</a>) con los de inestabilidad magnética quiral en magnetares de Dehman y Pons (<a href="https://arxiv.org/abs/2408.05281">arXiv:2408.05281</a>, <a href="https://arxiv.org/abs/2505.06196">arXiv:2505.06196</a>). Deduzco en SymPy desde la acción de Einstein–Maxwell–Axión que el gradiente radial de la pared de dominio del axión \\(a_0'(r) = f_a \\theta_0'(r)\\) actúa exactamente como un potencial químico quiral efectivo \\(\\mu_{{5,\\mathrm{{eff}}}}(r)\\) localizado en \\(r = r_{{\\mathrm{{wall}}}}\\), e integro el acoplo magneto-axiónico sobre el dipolo relativista de Wasserman–Shapiro en la estrella BSk26.</p>
    <div class="tally">
      <div><b>{checks_data["passed_checks"]}/{checks_data["total_checks"]}</b><span>comprobaciones OK</span></div>
      <div><b>{checks_data["negative_controls"]}</b><span>controles negativos superados</span></div>
      <div><b>0,87–1,24 cm</b><span>longitud de onda quiral \\(\\lambda_{{\\mathrm{{crit}}}}\\) en la pared</span></div>
      <div><b>0,022–0,031 eV</b><span>potencial quiral efectivo \\(\\mu_{{5,\\mathrm{{eff}}}}\\)</span></div>
      <div><b>{checks_data["elapsed_s"]:.1f} s</b><span>ejecución completa (8 estrellas × 4 campos)</span></div>
    </div>
  </header>

  <section class="stack" id="conexion">
    <h2>1. Qué conecta este cálculo entre las dos líneas del grupo</h2>
    <ol class="cites">
      <li>
        <span class="tag">Línea A · Condensados de axiones (2024–2026)</span>
        <p>Gómez-Bañón, Bartnick, Springmann, Pnigouras y Pons (<a href="https://arxiv.org/abs/2408.13247">arXiv:2408.13247</a>, <a href="https://arxiv.org/abs/2510.16245">2510.16245</a>, <a href="https://arxiv.org/abs/2610.37308">2610.37308</a>) estudian el condensado \\(\\theta_0(r) \\simeq \\pi\\) en el interior de una estrella de neutrones sin campo magnético (\\(\\mathbf{{B}} = 0\\)). Toda la dinámica nueva ocurre en la estrecha capa de transición \\(r \\approx r_{{\\mathrm{{wall}}}}\\) donde \\(\\theta_0(r)\\) cae de \\(\\pi\\) a \\(0\\).</p>
      </li>
      <li>
        <span class="tag">Línea B · Anomalía quiral y magnetares (2024–2025)</span>
        <p>Dehman y Pons (<a href="https://arxiv.org/abs/2408.05281">arXiv:2408.05281</a>, <a href="https://arxiv.org/abs/2505.06196">2505.06196</a>) estudian cómo una asimetría quiral \\(\\mu_5\\) en la corteza de una estrella de neutrones genera una corriente de efecto magnético quiral \\(\\mathbf{{J}}_{{\\mathrm{{CME}}}} = \\frac{{\\alpha_{{\\mathrm{{EM}}}}}}{{\\pi}}\\mu_5 \\mathbf{{B}}\\) que amplifica campos a gran escala con número de onda crítico \\(k_{{\\mathrm{{crit}}}} = \\frac{{2\\alpha_{{\\mathrm{{EM}}}}}}{{\\pi}}\\mu_5\\).</p>
      </li>
      <li>
        <span class="tag">Unificación · Acción de Einstein–Maxwell–Chern–Simons</span>
        <p>Al incluir el acoplo topológico estándar axión-fotón \\(-\\frac{{1}}{{4}} g_{{a\\gamma\\gamma}} a F_{{\\mu\\nu}}\\tilde{{F}}^{{\\mu\\nu}}\\) sobre la estrella condensada, la variación de Euler–Lagrange produce la corriente de Carroll–Field–Jackiw \\(J^\\nu_{{\\mathrm{{top}}}} = -g_{{a\\gamma\\gamma}}(\\partial_\\mu a)\\tilde{{F}}^{{\\mu\\nu}}\\). En la pared de dominio donde \\(\\partial_r a_0 = f_a \\theta_0'(r) \\neq 0\\), esta corriente tiene <b>exactamente la misma forma matemática que la corriente quiral de Dehman y Pons</b>, con un número de onda quiral localizado \\(k_{{\\mathrm{{chiral}}}}^{{\\mathrm{{ax}}}}(r) = |g_{{a\\gamma\\gamma}} f_a \\theta_0'(r)|\\) que no sufre lavado por colisiones de espín electrónicas.</p>
      </li>
    </ol>
  </section>

  <section class="stage wide" id="escena" aria-label="Escena interactiva Magnetar + Axión">
    <div class="stage-bar">
      <strong style="color:#fff;font-family:var(--f-display);font-size:1.08rem">Simulación interactiva: Dipolo relativista + pared de dominio quiral + conversión axión-fotón</strong>
      <button type="button" class="btn" id="playbtn" data-play>Pausar</button>
    </div>
    <div class="stage-canvas">
      <canvas id="magcanvas" width="960" height="400"></canvas>
    </div>
    <div class="stage-foot">
      <label class="slider">
        <span>Frecuencia del axión \\(\\nu_a\\): <b id="lbl-nu">6.0 kHz</b></span>
        <input type="range" id="sl-nu" min="0" max="7" step="1" value="3">
      </label>
      <label class="slider">
        <span>Campo polar \\(B_{{\\mathrm{{pole}}}}\\): <b id="lbl-b">10¹⁵ G</b></span>
        <input type="range" id="sl-b" min="0" max="3" step="1" value="3">
      </label>
      <div class="readout">
        <span>Pared \\(r_{{\\mathrm{{wall}}}}\\): <b id="ro-wall">10.43 km</b></span>
        <span>\\(\\mu_{{5,\\mathrm{{eff}}}}\\): <b id="ro-mu5">0.0298 eV</b></span>
        <span>\\(\\lambda_{{\\mathrm{{crit}}}}\\): <b id="ro-lam">0.89 cm</b></span>
        <span>\\(L_\\gamma\\): <b id="ro-lg">1.78×10⁴⁴ erg/s</b></span>
      </div>
    </div>
  </section>

  <section class="stack" id="derivacion">
    <h2>2. Derivación simbólica desde la acción en espacio-tiempo curvo</h2>
    <p>Partimos en <code>axion_ns/magnetar.py</code> de la densidad lagrangiana acoplada sobre la métrica esférica estática \\(ds^2 = e^{{2\\Phi_0}}dt^2 - e^{{2\\lambda_0}}dr^2 - r^2(d\\theta^2 + \\sin^2\\theta\\,d\\phi^2)\\):</p>
    \\[
      \\mathcal{{L}}_{{a+\\mathrm{{EM}}}} = \\sqrt{{-g}}\\left[\\frac{{1}}{{2}}g^{{\\mu\\nu}}\\partial_\\mu a\\,\\partial_\\nu a - V(a) - \\varepsilon(n,a) - \\frac{{1}}{{4}}F_{{\\mu\\nu}}F^{{\\mu\\nu}} - \\class{{newterm}}{{\\frac{{1}}{{4}}g_{{a\\gamma\\gamma}}\\,a\\,F_{{\\mu\\nu}}\\tilde{{F}}^{{\\mu\\nu}}}} - A_\\mu J^\\mu\\right],
    \\]
    <p>donde \\(\\tilde{{F}}^{{\\mu\\nu}} = \\frac{{1}}{{2\\sqrt{{-g}}}}\\varepsilon^{{\\mu\\nu\\alpha\\beta}}F_{{\\alpha\\beta}}\\). SymPy verifica primero las dos identidades topológicas exactas (\\(\\nabla_\\mu \\tilde{{F}}^{{\\mu\\nu}} \\equiv 0\\) y que \\(\\frac{{1}}{{4}}\\sqrt{{-g}}F_{{\\mu\\nu}}\\tilde{{F}}^{{\\mu\\nu}} = \\frac{{1}}{{2}}\\partial_\\mu(\\sqrt{{-g}}A_\\nu \\tilde{{F}}^{{\\mu\\nu}})\\)), gracias a las cuales las ecuaciones de Euler–Lagrange respecto a \\(A_\\nu\\) y \\(a\\) dependen únicamente de las derivadas primeras del axión:</p>
    \\[
      \\frac{{1}}{{\\sqrt{{-g}}}}\\partial_\\mu\\!\\left(\\sqrt{{-g}}\\,F^{{\\mu\\nu}}\\right) + \\class{{newterm}}{{g_{{a\\gamma\\gamma}}(\\partial_\\mu a)\\tilde{{F}}^{{\\mu\\nu}}}} = J^\\nu, \\qquad
      \\Box a + V'(a) + n\\frac{{\\partial m^*}}{{\\partial a}} + \\class{{newterm}}{{\\frac{{1}}{{4}}g_{{a\\gamma\\gamma}}F_{{\\mu\\nu}}\\tilde{{F}}^{{\\mu\\nu}}}} = 0.
    \\]

    <div class="ob">
      <header><span class="weight w-ok">Resultado 1 · Dipolo relativista y efecto Witten</span><h3>Ecuación de Wasserman–Shapiro (1983) y carga topológica en la pared</h3></header>
      <p>Para un campo magnético poloidal dipolar \\(A_\\phi(r,\\theta) = -a_\\phi(r)\\sin^2\\theta\\), la componente \\(\\nu=\\phi\\) reproduce exactamente la ecuación relativista de Wasserman y Shapiro:</p>
      \\[
        e^{{-2\\lambda_0}}a_\\phi'' + e^{{-2\\lambda_0}}(\\Phi_0' - \\lambda_0')a_\\phi' - \\frac{{2}}{{r^2}}a_\\phi = -4\\pi j_\\phi(r),
      \\]
      <p>mientras que la componente temporal \\(\\nu = t\\) muestra que el salto radial \\(a_0'(r) = f_a\\theta_0'(r)\\) en la pared de dominio induce una densidad de carga eléctrica topológica de Witten incluso en reposo:</p>
      \\[
        \\rho_{{\\mathrm{{top}}}}(r,\\theta) = e^{{\\Phi_0}}J^t_{{\\mathrm{{CS}}}} = \\class{{newterm}}{{-g_{{a\\gamma\\gamma}}f_a\\,e^{{-\\lambda_0}}\\theta_0'(r)\\,B^{{\\hat{{r}}}}_0(r,\\theta)}}.
      \\]
    </div>

    <div class="ob">
      <header><span class="weight w-main">Resultado 2 · Conexión directa con Dehman &amp; Pons (2024/2025)</span><h3>Potencial químico quiral efectivo \\(\\mu_{{5,\\mathrm{{eff}}}}(r)\\) en la pared del axión</h3></header>
      <p>Al perturbar el campo con una componente toroidal \\(\\delta A_\\theta = \\mathcal{{B}}(r)\\sin\\theta\\,e^{{-i\\omega t}}\\) y poloidal \\(\\delta A_\\phi = \\mathcal{{E}}(r)\\sin^2\\theta\\,e^{{-i\\omega t}}\\), el gradiente radial del condensado \\(a_0'(r) = f_a\\theta_0'(r)\\) introduce en la ecuación azimutal de Maxwell el término de corriente helicoidal:</p>
      \\[
        e^{{-2\\lambda_0}}\\mathcal{{E}}'' + e^{{-2\\lambda_0}}(\\Phi_0'-\\lambda_0')\\mathcal{{E}}' + \\left(\\omega^2 e^{{-2\\Phi_0}} - \\frac{{2}}{{r^2}}\\right)\\mathcal{{E}} + \\class{{newterm}}{{i\\omega\\,e^{{-(\\Phi_0+\\lambda_0)}}\\,g_{{a\\gamma\\gamma}}f_a\\theta_0'(r)\\,\\mathcal{{B}}(r)}} = 0.
      \\]
      <p>Comparando con la corriente magnética quiral \\(\\mathbf{{J}}_{{\\mathrm{{CME}}}} = \\frac{{\\alpha_{{\\mathrm{{EM}}}}}}{{\\pi}}\\mu_5\\mathbf{{B}}\\) de Dehman y Pons (<a href="https://arxiv.org/abs/2408.05281">arXiv:2408.05281</a>), identificamos de forma exacta el número de onda quiral topológico y el potencial químico quiral efectivo de la pared de dominio:</p>
      \\[
        \\class{{newterm}}{{k_{{\\mathrm{{chiral}}}}^{{\\mathrm{{ax}}}}(r) = g_{{a\\gamma\\gamma}}f_a\\theta_0'(r), \\qquad \\mu_{{5,\\mathrm{{eff}}}}(r) = \\frac{{\\pi}}{{\\alpha_{{\\mathrm{{EM}}}}}}g_{{a\\gamma\\gamma}}f_a\\theta_0'(r), \\qquad k_{{\\mathrm{{crit}}}}(r) = 2\\,g_{{a\\gamma\\gamma}}f_a\\theta_0'(r).}}
      \\]
    </div>

    <div class="ob">
      <header><span class="weight w-ok">Resultado 3 · Oscilaciones magneto-axiónicas acopladas</span><h3>Conversión coherente \\(\\delta\\theta(r) \\leftrightarrow \\mathcal{{E}}(r)\\) sobre el campo dipolar \\(B_0(r)\\)</h3></header>
      <p>Para una oscilación del condensado \\(a = a_0(r) + f_a\\delta\\theta(r)\\cos\\theta\\,e^{{-i\\omega t}}\\) sobre el campo magnético de fondo \\(B_0(r)\\), las ecuaciones de Klein–Gordon y Maxwell quedan acopladas por operadores adjuntos exactos:</p>
      \\[
        \\begin{{aligned}}
        e^{{-2\\lambda_0}}\\delta\\theta'' + e^{{-2\\lambda_0}}\\left(\\Phi_0'-\\lambda_0'+\\frac{{2}}{{r}}\\right)\\delta\\theta' + \\left(\\omega^2 e^{{-2\\Phi_0}} - m_{{\\mathrm{{eff}}}}^2 - \\frac{{2}}{{r^2}}\\right)\\delta\\theta - \\class{{newterm}}{{i\\omega\\frac{{g_{{a\\gamma\\gamma}}B_0(r)}}{{f_a r}}e^{{-\\Phi_0}}\\mathcal{{E}}}} &= 0, \\\\[4pt]
        e^{{-2\\lambda_0}}\\mathcal{{E}}'' + e^{{-2\\lambda_0}}(\\Phi_0'-\\lambda_0')\\mathcal{{E}}' + \\left(\\omega^2 e^{{-2\\Phi_0}} - \\frac{{2}}{{r^2}}\\right)\\mathcal{{E}} + \\class{{newterm}}{{i\\omega\\,g_{{a\\gamma\\gamma}}f_a B_0(r)\\,r\\,e^{{-\\Phi_0}}\\delta\\theta}} &= 0.
        \\end{{aligned}}
      \\]
    </div>
  </section>

  <section class="stack wide" id="graficas">
    <h2>3. Perfiles numéricos sobre la estrella BSk26 (compilados desde SymPy)</h2>
    <div class="charts">
      <article class="chart">
        <h4>A. Potencial quiral efectivo \\(\\mu_{{5,\\mathrm{{eff}}}}(r)\\) en la pared de dominio</h4>
        <p class="sub">Localización en la transición \\(\\theta_0: \\pi \\to 0\\) para \\(\\nu_a = 3, 6, 10, 15\\) kHz (\\(g_{{a\\gamma\\gamma}} = 10^{{-11}}\\,\\mathrm{{GeV}}^{{-1}}\\))</p>
        <svg class="plot" id="svg-chiral" viewBox="0 0 500 280"></svg>
      </article>
      <article class="chart">
        <h4>B. Campo dipolar relativista de Wasserman–Shapiro \\(B^{{\\hat{{r}}}}(r)/B_{{\\mathrm{{pole}}}}\\)</h4>
        <p class="sub">Amplificación por compacidad \\(2GM/R\\) y empalme interior-exterior exacto (error \\(5.1\\times 10^{{-9}}\\))</p>
        <svg class="plot" id="svg-dipole" viewBox="0 0 500 280"></svg>
      </article>
      <article class="chart">
        <h4>C. Luminosidad de conversión magneto-axiónica \\(L_\\gamma(B_{{\\mathrm{{pole}}}}, \\nu_a)\\)</h4>
        <p class="sub">Escala cuadrática exacta \\(L_\\gamma \\propto B_{{\\mathrm{{pole}}}}^2\\) para \\(B_{{\\mathrm{{pole}}}} \\in [10^{{13}}, 10^{{15}}]\\) G</p>
        <svg class="plot" id="svg-lum" viewBox="0 0 500 280"></svg>
      </article>
      <article class="chart">
        <h4>D. Desplazamiento de frecuencia magneto-axiónico \\(\\Delta\\nu(10^{{15}}\\,\\mathrm{{G}})\\)</h4>
        <p class="sub">Corrimiento de la cavidad de la pared de dominio en función de \\(\\nu_a\\)</p>
        <svg class="plot" id="svg-shift" viewBox="0 0 500 280"></svg>
      </article>
    </div>
  </section>

  <section class="stack wide" id="tabla">
    <h2>4. Tabla de resultados numéricos (estrella BSk26, \\(p_c = 100\\,\\mathrm{{MeV/fm^3}}\\))</h2>
    <p class="muted">Sin axión: \\(M = 1.536\\,M_\\odot\\), \\(R = 11.72\\,\\mathrm{{km}}\\), factor de amplificación dipolar relativista de Wasserman–Shapiro \\(f_{{\\mathrm{{WS}}}} = 1.3587\\).</p>
    <div class="tablewrap">
      <table class="steps">
        <thead>
          <tr>
            <th class="num">\\(\\nu_a\\) [kHz]</th>
            <th class="num">\\(M\\) [\\(M_\\odot\\)]</th>
            <th class="num">\\(R\\) [km]</th>
            <th class="num">\\(r_{{\\mathrm{{wall}}}}\\) [km]</th>
            <th class="num">Factor WS \\(f_{{\\mathrm{{WS}}}}\\)</th>
            <th class="num">\\(k_{{\\mathrm{{crit}}}}\\) [cm\\(^{{-1}}\\)]</th>
            <th class="num">\\(\\lambda_{{\\mathrm{{crit}}}}\\) [cm]</th>
            <th class="num">\\(\\mu_{{5,\\mathrm{{eff}}}}\\) [eV]</th>
            <th class="num">\\(\\Delta\\nu(10^{{15}}\\,\\mathrm{{G}})\\) [Hz]</th>
            <th class="num">\\(L_\\gamma(10^{{14}}\\,\\mathrm{{G}})\\) [erg/s]</th>
            <th class="num">\\(L_\\gamma(5\\cdot 10^{{14}}\\,\\mathrm{{G}})\\) [erg/s]</th>
            <th class="num">\\(L_\\gamma(10^{{15}}\\,\\mathrm{{G}})\\) [erg/s]</th>
          </tr>
        </thead>
        <tbody>
          {table_rows_str}
        </tbody>
      </table>
    </div>
  </section>

  <section class="stack wide" id="comprobaciones">
    <h2>5. Todas las comprobaciones automáticas ({checks_data["passed_checks"]}/{checks_data["total_checks"]} OK)</h2>
    <p class="muted">Incluye 4 controles negativos deliberados que mutan signos, factores métricos y multipolos para verificar que el comparador simbólico detecta cualquier desviación.</p>
    <ul class="checks">
      {checks_list_str}
    </ul>
  </section>

  <footer class="wide">
    <span>Proyecto <code>workshop_fisica</code> · Línea 3 (`axion_ns/magnetar.py`, `axion_ns/step13_magnetar.py`)</span>
    <span><a href="informe_axiones.html">← Volver al informe principal</a></span>
  </footer>
</main>

<script id="fig-data" type="application/json">{fig_json_inline}</script>
<script>
(() => {{
  const root = document.documentElement;
  const saved = localStorage.getItem('theme');
  if (saved) root.setAttribute('data-theme', saved);
  document.getElementById('themetoggle')?.addEventListener('click', () => {{
    const cur = root.getAttribute('data-theme') || (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    const next = cur === 'dark' ? 'light' : 'dark';
    root.setAttribute('data-theme', next);
    localStorage.setItem('theme', next);
  }});

  const DATA = JSON.parse(document.getElementById('fig-data').textContent);
  const rows = DATA.rows;
  const bKeys = ['1e+13', '1e+14', '5e+14', '1e+15'];
  const bLabels = ['10¹³ G', '10¹⁴ G', '5×10¹⁴ G', '10¹⁵ G'];

  // ---- Escena interactiva Canvas (estilo 3Blue1Brown)
  const cv = document.getElementById('magcanvas');
  const ctx = cv.getContext('2d');
  const slNu = document.getElementById('sl-nu');
  const slB = document.getElementById('sl-b');
  const playBtn = document.getElementById('playbtn');
  let playing = true, tAnim = 0;

  function updateReadout() {{
    const r = rows[+slNu.value] || rows[3];
    const bk = bKeys[+slB.value] || '1e+15';
    const binfo = r.magnetic_fields[bk];
    document.getElementById('lbl-nu').textContent = r.nu_a_khz.toFixed(1) + ' kHz';
    document.getElementById('lbl-b').textContent = bLabels[+slB.value];
    document.getElementById('ro-wall').textContent = r.r_wall_km.toFixed(2) + ' km';
    document.getElementById('ro-mu5').textContent = r.mu5_eff_max_ev.toFixed(4) + ' eV';
    document.getElementById('ro-lam').textContent = r.lambda_crit_cm.toFixed(2) + ' cm';
    document.getElementById('ro-lg').textContent = binfo.L_gamma_erg_s.toExponential(2) + ' erg/s';
  }}
  slNu.addEventListener('input', updateReadout);
  slB.addEventListener('input', updateReadout);
  playBtn.addEventListener('click', () => {{
    playing = !playing;
    playBtn.textContent = playing ? 'Pausar' : 'Reproducir';
  }});
  updateReadout();

  function drawStage() {{
    const W = cv.width, H = cv.height;
    ctx.clearRect(0, 0, W, H);
    const r = rows[+slNu.value] || rows[3];
    const bIdx = +slB.value;
    const cx = 290, cy = H / 2;
    const Rpx = 128 * (r.R_km / 11.72);
    const Rwall = Rpx * (r.r_wall_km / r.R_km);

    // Líneas de campo dipolar magnético relativista
    ctx.save();
    for (let L = 1; L <= 4; L++) {{
      const scale = Rpx * (0.85 + 0.42 * L);
      ctx.strokeStyle = `rgba(88, 196, 221, ${{0.18 + 0.08 * bIdx}})`;
      ctx.lineWidth = 1.2 + 0.35 * bIdx;
      for (const sgn of [-1, 1]) {{
        ctx.beginPath();
        for (let a = 0.08; a <= Math.PI - 0.08; a += 0.04) {{
          const rad = scale * Math.sin(a) * Math.sin(a);
          const x = cx + sgn * rad * Math.sin(a);
          const y = cy - rad * Math.cos(a);
          if (a === 0.08) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }}
        ctx.stroke();
      }}
    }}
    ctx.restore();

    // Núcleo condensado (theta_0 = pi)
    const grad = ctx.createRadialGradient(cx, cy, 4, cx, cy, Rpx);
    grad.addColorStop(0, 'rgba(124, 58, 237, 0.55)');
    grad.addColorStop(Rwall / Rpx, 'rgba(247, 217, 111, 0.75)');
    grad.addColorStop(1, 'rgba(30, 41, 59, 0.85)');
    ctx.fillStyle = grad;
    ctx.beginPath();
    ctx.arc(cx, cy, Rpx, 0, Math.PI * 2);
    ctx.fill();
    ctx.strokeStyle = '#64748b';
    ctx.lineWidth = 1.5;
    ctx.stroke();

    // Pared de dominio quiral vibrante (r = r_wall)
    ctx.strokeStyle = '#f7d96f';
    ctx.lineWidth = 2.8;
    ctx.beginPath();
    for (let a = 0; a <= Math.PI * 2 + 0.05; a += 0.05) {{
      const dr = 3.5 * Math.cos(a) * Math.sin(tAnim * 2.5);
      const x = cx + (Rwall + dr) * Math.cos(a);
      const y = cy + (Rwall + dr) * Math.sin(a);
      if (a === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }}
    ctx.stroke();

    // Pulsos de fotones convertidos en la pared de dominio
    const nPhot = 6 + 3 * bIdx;
    for (let k = 0; k < nPhot; k++) {{
      const ang = (k / nPhot) * Math.PI * 2 + 0.3;
      const prog = ((tAnim * 0.7 + k / nPhot) % 1);
      const rad = Rwall + prog * 130;
      const px = cx + rad * Math.cos(ang);
      const py = cy + rad * Math.sin(ang);
      ctx.fillStyle = `rgba(247, 217, 111, ${{1 - prog}})`;
      ctx.beginPath();
      ctx.arc(px, py, 2.5 + 0.6 * bIdx, 0, Math.PI * 2);
      ctx.fill();
    }}

    // Panel derecho: perfil radial de theta_0(r) y pico quiral mu_{{5,eff}}(r)
    const px0 = 560, py0 = 60, pw = 350, ph = 280;
    ctx.strokeStyle = '#2c3850';
    ctx.strokeRect(px0, py0, pw, ph);
    ctx.fillStyle = '#aab6c9';
    ctx.font = '12px "IBM Plex Sans", sans-serif';
    ctx.fillText('Perfil en la pared de dominio: θ₀(r)/π (azul) y corriente quiral J_top (oro)', px0, py0 - 12);

    const wFrac = r.r_wall_km / (1.3 * r.R_km);
    // Curva theta_0 / pi
    ctx.strokeStyle = '#58c4dd';
    ctx.lineWidth = 2.2;
    ctx.beginPath();
    for (let i = 0; i <= 120; i++) {{
      const u = i / 120;
      const thVal = 1 / (1 + Math.exp((u - wFrac) / 0.025));
      const x = px0 + u * pw;
      const y = py0 + ph - thVal * (ph - 35) - 15;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }}
    ctx.stroke();

    // Pico de corriente quiral |theta_0'(r)|
    ctx.strokeStyle = '#f7d96f';
    ctx.lineWidth = 2.4;
    ctx.beginPath();
    for (let i = 0; i <= 120; i++) {{
      const u = i / 120;
      const ex = Math.exp((u - wFrac) / 0.025);
      const dVal = (4 * ex / ((1 + ex) * (1 + ex))) * (0.55 + 0.14 * bIdx) * (0.85 + 0.15 * Math.sin(tAnim * 3));
      const x = px0 + u * pw;
      const y = py0 + ph - dVal * (ph - 35) - 15;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }}
    ctx.stroke();

    if (playing) tAnim += 0.03;
    requestAnimationFrame(drawStage);
  }}
  requestAnimationFrame(drawStage);

  // ---- Gráficas SVG
  const NS = 'http://www.w3.org/2000/svg';
  function el(tag, attrs, text) {{
    const n = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs || {{}})) n.setAttribute(k, v);
    if (text !== undefined) n.textContent = text;
    return n;
  }}

  // Gráfica A: Perfiles quirales mu_{{5,eff}}(r)
  (function plotChiral() {{
    const svg = document.getElementById('svg-chiral');
    const pad = {{ l: 52, r: 18, t: 18, b: 36 }}, W = 500, H = 280;
    const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
    const xMin = 6.0, xMax = 13.5, yMin = 0, yMax = 0.034;
    const sx = x => pad.l + ((x - xMin) / (xMax - xMin)) * pw;
    const sy = y => pad.t + ph - ((y - yMin) / (yMax - yMin)) * ph;

    for (let yv = 0; yv <= 0.03; yv += 0.01) {{
      svg.appendChild(el('line', {{ x1: pad.l, x2: W - pad.r, y1: sy(yv), y2: sy(yv), stroke: 'var(--rule)' }}));
      svg.appendChild(el('text', {{ x: pad.l - 6, y: sy(yv) + 4, 'text-anchor': 'end', fill: 'var(--muted)', 'font-size': 10 }}, yv.toFixed(2)));
    }}
    for (let xv = 6; xv <= 13; xv += 2) {{
      svg.appendChild(el('text', {{ x: sx(xv), y: H - 10, 'text-anchor': 'middle', fill: 'var(--muted)', 'font-size': 10 }}, xv + ' km'));
    }}
    const colors = {{ '3.0': 'var(--c1)', '6.0': 'var(--c2)', '10.0': 'var(--c3)', '15.0': 'var(--c4)' }};
    for (const [k, prof] of Object.entries(DATA.chiral_profiles)) {{
      const pts = prof.r_km.map((rk, i) => [rk, prof.mu5_eff_ev[i]]).filter(p => p[0] >= xMin && p[0] <= xMax);
      const d = pts.map((p, i) => (i ? 'L' : 'M') + sx(p[0]).toFixed(1) + ',' + sy(p[1]).toFixed(1)).join(' ');
      svg.appendChild(el('path', {{ d, fill: 'none', stroke: colors[k] || 'var(--c1)', 'stroke-width': 2.2 }}));
    }}
  }})();

  // Gráfica B: Dipolo relativista B_r(r) / B_pole
  (function plotDipole() {{
    const svg = document.getElementById('svg-dipole');
    const pad = {{ l: 52, r: 18, t: 18, b: 36 }}, W = 500, H = 280;
    const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
    const xMin = 0, xMax = 22, yMin = 0, yMax = 5.0;
    const sx = x => pad.l + ((x - xMin) / (xMax - xMin)) * pw;
    const sy = y => pad.t + ph - ((y - yMin) / (yMax - yMin)) * ph;

    for (let yv = 0; yv <= 5; yv += 1) {{
      svg.appendChild(el('line', {{ x1: pad.l, x2: W - pad.r, y1: sy(yv), y2: sy(yv), stroke: 'var(--rule)' }}));
      svg.appendChild(el('text', {{ x: pad.l - 6, y: sy(yv) + 4, 'text-anchor': 'end', fill: 'var(--muted)', 'font-size': 10 }}, yv.toFixed(0)));
    }}
    for (let xv = 0; xv <= 20; xv += 5) {{
      svg.appendChild(el('text', {{ x: sx(xv), y: H - 10, 'text-anchor': 'middle', fill: 'var(--muted)', 'font-size': 10 }}, xv + ' km'));
    }}
    const keys = [['no_axion', 'var(--muted)', '4,3'], ['3.0', 'var(--c1)', ''], ['6.0', 'var(--c2)', ''], ['15.0', 'var(--c3)', '']];
    for (const [k, col, dash] of keys) {{
      const prof = DATA.dipole_profiles[k];
      if (!prof) continue;
      const d = prof.r_km.map((rk, i) => (i ? 'L' : 'M') + sx(rk).toFixed(1) + ',' + sy(prof.B_r_over_Bpole[i]).toFixed(1)).join(' ');
      const attrs = {{ d, fill: 'none', stroke: col, 'stroke-width': 2.1 }};
      if (dash) attrs['stroke-dasharray'] = dash;
      svg.appendChild(el('path', attrs));
    }}
  }})();

  // Gráfica C: Luminosidad L_gamma vs nu_a
  (function plotLum() {{
    const svg = document.getElementById('svg-lum');
    const pad = {{ l: 56, r: 18, t: 18, b: 36 }}, W = 500, H = 280;
    const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
    const xMin = 2.5, xMax = 15.5, yMin = 39.5, yMax = 45.0;
    const sx = x => pad.l + ((x - xMin) / (xMax - xMin)) * pw;
    const sy = y => pad.t + ph - ((y - yMin) / (yMax - yMin)) * ph;

    for (let yv = 40; yv <= 45; yv += 1) {{
      svg.appendChild(el('line', {{ x1: pad.l, x2: W - pad.r, y1: sy(yv), y2: sy(yv), stroke: 'var(--rule)' }}));
      svg.appendChild(el('text', {{ x: pad.l - 6, y: sy(yv) + 4, 'text-anchor': 'end', fill: 'var(--muted)', 'font-size': 10 }}, '10^' + yv));
    }}
    for (const xv of [3, 6, 9, 12, 15]) {{
      svg.appendChild(el('text', {{ x: sx(xv), y: H - 10, 'text-anchor': 'middle', fill: 'var(--muted)', 'font-size': 10 }}, xv + ' kHz'));
    }}
    const bCols = [['1e+13', 'var(--c4)'], ['1e+14', 'var(--c1)'], ['5e+14', 'var(--c3)'], ['1e+15', 'var(--c2)']];
    for (const [bk, col] of bCols) {{
      const d = rows.map((r, i) => {{
        const lg = Math.log10(r.magnetic_fields[bk].L_gamma_erg_s);
        return (i ? 'L' : 'M') + sx(r.nu_a_khz).toFixed(1) + ',' + sy(lg).toFixed(1);
      }}).join(' ');
      svg.appendChild(el('path', {{ d, fill: 'none', stroke: col, 'stroke-width': 2.2 }}));
      for (const r of rows) {{
        const lg = Math.log10(r.magnetic_fields[bk].L_gamma_erg_s);
        svg.appendChild(el('circle', {{ cx: sx(r.nu_a_khz), cy: sy(lg), r: 3.5, fill: col }}));
      }}
    }}
  }})();

  // Gráfica D: Desplazamiento de frecuencia Delta nu (Hz) a 10^15 G
  (function plotShift() {{
    const svg = document.getElementById('svg-shift');
    const pad = {{ l: 56, r: 18, t: 18, b: 36 }}, W = 500, H = 280;
    const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
    const xMin = 2.5, xMax = 15.5, yMin = -200, yMax = 1250;
    const sx = x => pad.l + ((x - xMin) / (xMax - xMin)) * pw;
    const sy = y => pad.t + ph - ((y - yMin) / (yMax - yMin)) * ph;

    for (const yv of [0, 400, 800, 1200]) {{
      svg.appendChild(el('line', {{ x1: pad.l, x2: W - pad.r, y1: sy(yv), y2: sy(yv), stroke: 'var(--rule)' }}));
      svg.appendChild(el('text', {{ x: pad.l - 6, y: sy(yv) + 4, 'text-anchor': 'end', fill: 'var(--muted)', 'font-size': 10 }}, yv + ' Hz'));
    }}
    for (const xv of [3, 6, 9, 12, 15]) {{
      svg.appendChild(el('text', {{ x: sx(xv), y: H - 10, 'text-anchor': 'middle', fill: 'var(--muted)', 'font-size': 10 }}, xv + ' kHz'));
    }}
    const d = rows.map((r, i) => (i ? 'L' : 'M') + sx(r.nu_a_khz).toFixed(1) + ',' + sy(r.magnetic_fields['1e+15'].delta_nu_hz).toFixed(1)).join(' ');
    svg.appendChild(el('path', {{ d, fill: 'none', stroke: 'var(--c2)', 'stroke-width': 2.4 }}));
    for (const r of rows) {{
      svg.appendChild(el('circle', {{ cx: sx(r.nu_a_khz), cy: sy(r.magnetic_fields['1e+15'].delta_nu_hz), r: 4, fill: 'var(--c2)' }}));
    }}
  }})();
}})();
</script>
</body>
</html>
"""
    (REP / 'informe_magnetar.html').write_text(page)
    (RES / 'informe_magnetar.html').write_text(page)
    print(RES / 'informe_magnetar.html', f'{len(page) / 1e3:.0f} kB')


if __name__ == '__main__':
    build()
