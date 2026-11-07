"""Genera report/informe_radial.html (y la copia en results/) desde results/radial_branches.json.

Regla (AGENTS.md §4): ninguna cifra se teclea a mano; todas salen del JSON que escribe
axion_ns/step14_radial_branches.py.
"""
import html
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
D = json.loads((ROOT / 'results' / 'radial_branches.json').read_text())


def fnum(x, nd=3):
    return f'{x:.{nd}f}'.replace('-', '−')


# ------------------------------------------------------------------ SVG helpers
class Plot:
    def __init__(self, xr, yr, w=720, h=380, m=(16, 20, 46, 62), xlog=False, ylog=False):
        self.xr, self.yr, self.w, self.h = xr, yr, w, h
        self.mt, self.mr, self.mb, self.ml = m
        self.xlog, self.ylog = xlog, ylog
        self.parts = []

    def X(self, x):
        a, b = self.xr
        if self.xlog:
            x, a, b = math.log10(x), math.log10(a), math.log10(b)
        return self.ml + (x - a) / (b - a) * (self.w - self.ml - self.mr)

    def Y(self, y):
        a, b = self.yr
        if self.ylog:
            y, a, b = math.log10(y), math.log10(a), math.log10(b)
        return self.h - self.mb - (y - a) / (b - a) * (self.h - self.mt - self.mb)

    def axes(self, xt, yt, xlabel, ylabel, xfmt=str, yfmt=str):
        p = self.parts
        for v in yt:
            y = self.Y(v)
            p.append(f'<line class="grid" x1="{self.ml}" x2="{self.w - self.mr}" y1="{y:.1f}" y2="{y:.1f}"/>')
            p.append(f'<text class="tick" x="{self.ml - 8}" y="{y + 4:.1f}" text-anchor="end">{yfmt(v)}</text>')
        for v in xt:
            x = self.X(v)
            p.append(f'<line class="grid" x1="{x:.1f}" x2="{x:.1f}" y1="{self.mt}" y2="{self.h - self.mb}"/>')
            p.append(f'<text class="tick" x="{x:.1f}" y="{self.h - self.mb + 18}" text-anchor="middle">{xfmt(v)}</text>')
        p.append(f'<rect class="frame" x="{self.ml}" y="{self.mt}" width="{self.w - self.ml - self.mr}" height="{self.h - self.mt - self.mb}"/>')
        p.append(f'<text class="axlabel" x="{(self.ml + self.w - self.mr) / 2:.1f}" y="{self.h - 8}" text-anchor="middle">{xlabel}</text>')
        p.append(f'<text class="axlabel" transform="translate(16 {(self.mt + self.h - self.mb) / 2:.1f}) rotate(-90)" text-anchor="middle">{ylabel}</text>')

    def hline(self, y, cls, label=None, anchor='end'):
        yy = self.Y(y)
        self.parts.append(f'<line class="{cls}" x1="{self.ml}" x2="{self.w - self.mr}" y1="{yy:.1f}" y2="{yy:.1f}"/>')
        if label:
            x = self.w - self.mr - 6 if anchor == 'end' else self.ml + 6
            self.parts.append(f'<text class="note" x="{x}" y="{yy - 6:.1f}" text-anchor="{anchor}">{label}</text>')

    def vline(self, x, cls):
        xx = self.X(x)
        self.parts.append(f'<line class="{cls}" x1="{xx:.1f}" x2="{xx:.1f}" y1="{self.mt}" y2="{self.h - self.mb}"/>')

    def path(self, pts, cls):
        d = ' '.join(f'{"M" if i == 0 else "L"}{self.X(x):.1f},{self.Y(y):.1f}' for i, (x, y) in enumerate(pts))
        self.parts.append(f'<path class="{cls}" d="{d}"/>')

    def dots(self, pts, cls, tip, r=4):
        for x, y in pts:
            self.parts.append(f'<circle class="{cls}" cx="{self.X(x):.1f}" cy="{self.Y(y):.1f}" r="{r}"><title>{tip(x, y)}</title></circle>')

    def text(self, x, y, s, cls='note', anchor='start', dx=0, dy=0):
        self.parts.append(f'<text class="{cls}" x="{self.X(x) + dx:.1f}" y="{self.Y(y) + dy:.1f}" text-anchor="{anchor}">{s}</text>')

    def svg(self, label):
        return (f'<svg viewBox="0 0 {self.w} {self.h}" role="img" aria-label="{html.escape(label)}">'
                + ''.join(self.parts) + '</svg>')


# ------------------------------------------------------------------ Figura 1: plano k
tr = D['threshold_track']
after = D['axis_after']
x_coll = tr[-1][0]
x_prev = tr[-2][0]
pair = [r for r in tr if abs(r[1]) > 1e-9]
b_lo, b_hi = -pair[0][2], -pair[-1][2]
bound86 = D['axis_x0_0p86']['bound'][0]
virt86 = D['axis_x0_0p86']['virtual'][0]
w86 = math.sqrt(1 - bound86**2)

p1 = Plot((-0.42, 0.42), (-0.2, 0.12), h=360)
p1.axes([-0.4, -0.2, 0, 0.2, 0.4], [-0.2, -0.1, 0, 0.1], 'Re k', 'Im k', lambda v: fnum(v, 1), lambda v: fnum(v, 1))
p1.hline(0, 'zero')
p1.vline(0, 'zero')
p1.path([(r[1], r[2]) for r in pair], 's1-line')
p1.path([(-r[1], r[2]) for r in pair], 's1-line')
tipq = lambda x, y: f'cuasi-normal: k = {fnum(x, 4)} {"+" if y >= 0 else "−"} {fnum(abs(y), 4)} i'
p1.dots([(r[1], r[2]) for r in pair[::4]] + [(-r[1], r[2]) for r in pair[::4]], 's1-dot', tipq, r=3.5)
bd = [(0, k) for a in after for k in a['kappa'] if k > 0]
vd = [(0, k) for a in after for k in a['kappa'] if k < 0]
p1.dots(bd, 's3-dot', lambda x, y: f'ligado: k = +{fnum(y, 4)} i', r=4)
p1.dots(vd, 's2-dot', lambda x, y: f'virtual: k = −{fnum(-y, 4)} i', r=4)
p1.text(pair[0][1], pair[0][2], f'x₀ = {fnum(pair[0][0], 2)}', dy=-10, anchor='middle')
p1.text(0, max(k for _, k in bd), 'ligado (Im k > 0)', dx=10, dy=4)
p1.text(0, min(k for _, k in vd), 'virtual (Im k < 0)', dx=10, dy=4)
p1.text(-0.40, 0.085, 'evanescente fuera: estados ligados', cls='zone')
p1.text(-0.40, -0.185, 'onda saliente o creciente', cls='zone')
fig1 = p1.svg('Trayectoria de las raíces en el plano k al variar x0')

# ------------------------------------------------------------------ Figura 2: ramas w^2(x0)
p2 = Plot((0.80, 1.30), (-1.2, 1.1), h=380)
p2.axes([0.8, 0.9, 1.0, 1.1, 1.2, 1.3], [-1, -0.5, 0, 0.5, 1], 'x₀ = m_a R', 'w² = (ω/m_a)²',
        lambda v: fnum(v, 1), lambda v: fnum(v, 1))
p2.hline(1, 'thr', 'umbral ω = m_a')
p2.hline(0, 'unst', 'ω² = 0: por debajo, modo que crece', anchor='start')
pts0 = [(x0, w2) for x0, ws in D['branches']['0.0'] for w2 in ws if -1.2 < w2 < 1.1]
pts1 = [(x0, w2) for x0, ws in D['branches']['0.001'] for w2 in ws if -1.2 < w2 < 1.1]
p2.dots(pts0, 'ref-dot', lambda x, y: f'λ² = 0, x₀ = {fnum(x, 2)}: w² = {fnum(y, 4)}', r=3.5)
p2.dots(pts1, 's1-dot', lambda x, y: f'λ² = 10⁻³, x₀ = {fnum(x, 2)}: w² = {fnum(y, 4)}', r=3.5)
deep = [w2 for _, ws in D['branches']['0.001'] for w2 in ws if w2 < -1.2]
gap_x0 = D['gaps'][-1]['x0']
p2.text(gap_x0 + 0.012, 1.0, '← cruce evitado', anchor='start', dy=-8)
fig2 = p2.svg('Estados ligados frente a x0, con y sin acoplo')

# ------------------------------------------------------------------ Figura 3: amortiguamiento cerca del umbral
dp = D['damping_points']
s_fit = D['damping_exponent']
xmin, xmax = min(a for a, _ in dp), max(a for a, _ in dp)
ymin, ymax = min(b for _, b in dp), max(b for _, b in dp)
p3 = Plot((10 ** math.floor(math.log10(xmin)), 10 ** math.ceil(math.log10(xmax))),
          (10 ** math.floor(math.log10(ymin)), 10 ** math.ceil(math.log10(ymax))), h=330, xlog=True, ylog=True)
xt = [10 ** e for e in range(int(math.log10(p3.xr[0])), int(math.log10(p3.xr[1])) + 1)]
yt = [10 ** e for e in range(int(math.log10(p3.yr[0])), int(math.log10(p3.yr[1])) + 1)]
sup = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
p3.axes(xt, yt, 'Re w − 1', '|Im w|', lambda v: '10' + str(round(math.log10(v))).translate(sup),
        lambda v: '10' + str(round(math.log10(v))).translate(sup))
xm, ym = dp[len(dp) // 2]
ref = [(x, ym * (x / xm) ** 0.5) for x in (xmin, xmax)]
p3.path(ref, 'ref-line')
p3.dots(dp[::2], 's1-dot', lambda x, y: f'Re w − 1 = {x:.2e}, |Im w| = {y:.2e}', r=3.5)
p3.text(xmax, ref[1][1], 'pendiente ½', anchor='end', dy=-10)
fig3 = p3.svg('Amortiguamiento frente a la distancia al umbral, escala log-log')

# ------------------------------------------------------------------ tablas
gap_rows = ''.join(
    f'<tr><td>10<sup>{round(math.log10(g["lam2"]))}</sup></td><td>{g["gap_w2"]:.5f}</td>'
    f'<td>{g["gap_w2"] / math.sqrt(g["lam2"]):.3f}</td><td>{g["x0"]:.4f}</td></tr>' for g in D['gaps'])
ratios = [g['gap_w2'] / math.sqrt(g['lam2']) for g in D['gaps']]
checks = D['checks']
n_ok = sum(c['ok'] for c in checks)
check_rows = ''.join(
    f'<li class="{"ok" if c["ok"] else "bad"}"><span class="pill">{"OK" if c["ok"] else "FALLA"}</span>'
    f'<span>{html.escape(c["name"])}{(" <em>(" + html.escape(c["note"]) + ")</em>") if c["note"] else ""}</span></li>'
    for c in checks)

# estimación a priori para H11 (no es un resultado: se calcula aquí para no teclearla)
GM_SUN_KM, C_KMS = 1.476625, 299792.458
alpha_g = 1.5 * GM_SUN_KM * 2 * math.pi * 3000 / C_KMS
wf = D['fluid_x0_0p6']
xc, xc_an = D['x0_unstable'], D['x0_unstable_analytic']
w2u = D['unstable_x0_0p6_w2']
grow = math.sqrt(-w2u)
alpha = D['params']['alpha']

# ------------------------------------------------------------------ estrella realista (paso 15)
SP = ROOT / 'results' / 'radial_star.json'
S = json.loads(SP.read_text()) if SP.exists() else None
if S:
    s6 = S['nu6']
    # Figura 4: espectro a nu_a = 6 kHz en dos filas (frecuencia real / amortiguados)
    p4 = Plot((3.0, 9.0), (-1.6, 1.6), h=230, m=(16, 20, 46, 20))
    p4.parts.append(f'<rect class="frame" x="{p4.ml}" y="{p4.mt}" width="{p4.w - p4.ml - p4.mr}" height="{p4.h - p4.mt - p4.mb}"/>')
    for v in (3, 4, 5, 6, 7, 8, 9):
        p4.parts.append(f'<text class="tick" x="{p4.X(v):.1f}" y="{p4.h - p4.mb + 18}" text-anchor="middle">{v}</text>')
    p4.parts.append(f'<text class="axlabel" x="{(p4.ml + p4.w - p4.mr) / 2:.1f}" y="{p4.h - 8}" text-anchor="middle">Re ν (kHz)</text>')
    p4.hline(0.8, 'grid'); p4.hline(-0.8, 'grid')
    p4.vline(6.0, 'thr')
    p4.text(6.0, 1.35, 'ν_a = 6 kHz', anchor='start', dx=6)
    p4.text(8.95, 0.8, 'frecuencia real (atrapados)', cls='zone', anchor='end', dy=-10)
    p4.text(3.05, -0.8, 'amortiguados (Im ν &lt; 0)', cls='zone', dy=-14)
    p4.dots([(s6['F'], 0.8)], 's1-dot', lambda x, y: f'modo F: {x:.4f} kHz, real', r=6)
    lad = [l['nu'] for l in s6['ladder']]
    p4.dots([(x, 0.8) for x in lad], 's2-dot', lambda x, y: f'átomo gravitatorio: {x:.5f} kHz', r=4)
    p4.dots([(s6['H1'][0], -0.8)], 's1-ring', lambda x, y: f'H1: {s6["H1"][0]:.4f} {s6["H1"][1] * 1e3:+.1f}i Hz', r=6)
    p4.dots([(s6['axion_res'][0], -0.8)], 's2-ring', lambda x, y: f'resonancia axiónica: {s6["axion_res"][0]:.3f} {s6["axion_res"][1] * 1e3:+.0f}i Hz', r=6)
    p4.text(s6['F'], 0.8, f'F {fnum(s6["F"], 3)}', anchor='middle', dy=24)
    p4.text(min(lad), 0.8, 'átomo n ≥ 2', anchor='end', dx=-10, dy=4)
    p4.text(s6['H1'][0], -0.8, f'H1, Im ν = {fnum(s6["H1"][1] * 1e3, 1)} Hz', anchor='middle', dy=22)
    p4.text(s6['axion_res'][0], -0.8, f'axión, Im ν = {fnum(s6["axion_res"][1] * 1e3, 0)} Hz', anchor='middle', dy=22)
    fig4 = p4.svg('Espectro radial a nu_a = 6 kHz')
    # Figura 5: escalera hidrogenoide
    ag = s6['alpha_g']
    p5 = Plot((1.5, 6.5), (1e-4, 1e-2), h=300, ylog=True)
    p5.axes([2, 3, 4, 5, 6], [1e-4, 1e-3, 1e-2], 'n', '1 − ν/ν_a', lambda v: str(v),
            lambda v: '10' + str(round(math.log10(v))).translate(sup))
    p5.path([(n / 10, ag**2 / (2 * (n / 10)**2)) for n in range(15, 66)], 'ref-line')
    p5.dots([(l['n'], 1 - l['nu'] / 6.0) for l in s6['ladder']], 's2-dot',
            lambda x, y: f'n = {round(x)}: 1 − ν/ν_a = {y:.3e}', r=5)
    p5.text(2.6, ag**2 / (2 * 2.6**2), f'α_g²/(2n²), α_g = {fnum(ag, 4)}', dx=8, dy=-8)
    fig5 = p5.svg('Niveles del átomo gravitatorio frente a la ley hidrogenoide')
    lad_rows = ''.join(f'<tr><td>{l["n"]}</td><td>{l["nu"]:.6f}</td><td>{l["n_eff"]:.3f}</td></tr>' for l in s6['ladder'])
    t = S['turning']
    n35, n34 = S['nu35'], S['nu34']
    st_checks = S.get('checks', [])
    st_rows = ''.join(
        f'<li class="{"ok" if c["ok"] else "bad"}"><span class="pill">{"OK" if c["ok"] else "FALLA"}</span>'
        f'<span>{html.escape(c["name"])}{(" <em>(" + html.escape(c["note"]) + ")</em>") if c["note"] else ""}</span></li>'
        for c in st_checks)
    star_html = f"""
<section>
  <h2>6. La estrella realista: BSk26 con condensado, en relatividad general completa</h2>
  <p>Las ecuaciones radiales salen de la misma linealización de la acción que el paso 2 compara con la Ec. (3.3) del Paper I. Se convierten en el sistema <code>(ξ, Δp, δθ, δθ')</code> sin transcribir nada a mano (<code>axion_ns/radial.py</code>):</p>
  <ul>
    <li><code>δλ</code> y <code>δλ'</code> salen de <code>G_tr</code>, <code>δΦ'</code> de <code>G_rr</code> y <code>ξ'</code> de la conservación bariónica.</li>
    <li>El condensado no se anula en la superficie, así que fuera de la estrella se sigue integrando el axión acoplado a la métrica hasta que la cola es despreciable.</li>
    <li>Más allá, la onda exterior se integra en un rayo del plano r complejo donde la onda saliente decae. En la recta real, la integración hacia dentro es inestable para modos muy amortiguados: el método anterior (<code>axstar.exterior_L</code>) se desvía un {fnum(S['exterior']['old_vs_new_strong'] * 100, 0)} % en la resonancia axiónica, aunque solo {S['exterior']['old_vs_new_weak']:.0e} en un modo poco amortiguado, así que los amortiguamientos de pocos Hz de los pasos 9–10 siguen siendo válidos.</li>
  </ul>
  <p><b>Validación sin axión.</b> El cuadrado de la frecuencia del modo F se anula en <code>p_c = {fnum(t['p0_wF2_zero'], 1)}</code> MeV/fm³ y la masa máxima ({fnum(t['Mmax_msun'], 4)} M☉) está en <code>p_c = {fnum(t['p0_Mmax'], 1)}</code>. Coinciden, como exige el criterio del punto de retorno.</p>
  <figure>
    <div class="legend"><span><i class="sw s1"></i>modo fluido</span><span><i class="sw s2"></i>modo axiónico</span><span><i class="sw ref"></i>anillo: amortiguado</span></div>
    {fig4}
    <figcaption>Espectro radial a ν_a = 6 kHz (M = {fnum(s6['M_msun'], 4)} M☉, R = {fnum(s6['R_km'], 3)} km). El carácter se mide con la fracción de energía en el fluido: F {fnum(s6['F_fluid'], 2)}, H1 {fnum(s6['H1_fluid'], 2)}, resonancia axiónica {fnum(s6['axion_res_fluid'], 2)}.</figcaption>
  </figure>
  <p>La estructura es la del modelo simplificado, con un añadido de la gravedad:</p>
  <ul>
    <li>Los modos fluidos por debajo del umbral tienen frecuencia real (F) y los de encima se amortiguan poco (H1).</li>
    <li>La rama axiónica da resonancias muy amortiguadas.</li>
    <li>Justo por debajo del umbral, el potencial gravitatorio exterior atrapa al axión en una <b>escalera hidrogenoide</b>: un «átomo gravitatorio».</li>
  </ul>
  <figure>
    {fig5}
    <figcaption>Energía de ligadura de los niveles axiónicos frente a la ley hidrogenoide <code>α_g²/(2n²)</code>, con <code>α_g = GMm_a</code> (sin ajuste). El nivel n = 1, cuya órbita se solapa con la estrella, se desvía; se trata en el cuaderno.</figcaption>
  </figure>
  <div class="tbl"><table><tr><th>n</th><th>ν (kHz)</th><th>n_eff = α_g/√(2(1 − ν/ν_a))</th></tr>{lad_rows}</table></div>
  <h3>Cómo cruza el modo F el umbral</h3>
  <p>Al bajar ν_a, el umbral se acerca a F desde arriba, y el nivel n = 1 del átomo gravitatorio baja con él:</p>
  <ul>
    <li>En <b>ν_a = 3,5 kHz</b> hay dos estados ligados ({fnum(n35['bound'][0], 4)} y {fnum(n35['bound'][1], 4)} kHz) con carácter mezclado ({fnum(n35['fluid'][0], 2)} y {fnum(n35['fluid'][1], 2)} de fluido): un cruce evitado entre F y el átomo.</li>
    <li>En <b>ν_a = 3,4 kHz</b>, F ya está por encima del umbral, en {fnum(n34['F'][0], 4)} kHz, con un amortiguamiento de solo {fnum(-n34['F'][1] * 1e3, 2)} Hz (fracción fluida {fnum(n34['fluid'], 2)}).</li>
  </ul>
  <p>La transición cuasi-normal → ligado de la estrella realista no es un cruce limpio como en el modelo simplificado: el modo fluido atraviesa una escalera densa de estados del axión que se acumula en el umbral.</p>
  <div class="callout">
    <h3>Abierto</h3>
    <p>No hemos podido decidir si el amortiguamiento tiende a cero o a un valor finito justo encima del umbral (factor de Sommerfeld del potencial 1/r). Los dos puntos disponibles apuntan a que cae deprisa, pero hace falta un barrido denso en ν_a.</p>
  </div>
  <ul class="checks">{st_rows}</ul>
</section>"""
    star_summary = f"""
    <li><b>Estrella realista</b> (BSk26, relatividad general completa, ν_a = 6 kHz): el modo F por debajo del umbral tiene frecuencia real, H1 por encima se amortigua {fnum(-s6['H1'][1] * 1e3, 1)} Hz y la rama axiónica da resonancias muy amortiguadas. Justo bajo el umbral aparece una escalera de estados ligados del axión que sigue la ley hidrogenoide, con <code>n_eff</code> = {', '.join(f"{l['n_eff']:.2f}" for l in s6['ladder'])}: un «átomo gravitatorio».</li>
    <li>En la estrella realista, el modo F cruza el umbral atravesando esa escalera, con un cruce evitado con el nivel n = 1, y sale como cuasi-normal con un amortiguamiento muy pequeño ({fnum(-n34['F'][1] * 1e3, 2)} Hz a ν_a = 3,4 kHz).</li>"""
else:
    star_html = None
    star_summary = ''

PENDING = f'''<section>
  <h2>6. Siguiente paso: la estrella realista <span class="status pend">· pendiente</span></h2>
  <p>Integraremos los modos radiales en relatividad general completa sobre la estrella BSk26 con condensado, a partir de las Ec. (3.3) del Paper I verificadas en el paso 2. Predicciones anotadas antes de calcular:</p>
  <ol>
    <li><b>Validación sin axión:</b> <code>ω_F² → 0</code> en la configuración de masa máxima (criterio de punto de retorno).</li>
    <li>Con axión, <code>Im ν = 0</code> para los modos con <code>ν &lt; ν_a</code> e <code>Im ν &lt; 0</code> por encima.</li>
    <li><b>Especulativa:</b> la gravedad exterior <code>1/r</code> crea una escalera hidrogenoide de estados axiónicos ligados que se acumula en el umbral, <code>ω_n ≈ m_a[1 − α_g²/(2n²)]</code>. El acoplo gravitatorio es <code>α_g = GMm_a ≈ {fnum(alpha_g, 2)}</code> para 1,5 M☉ y 3 kHz (estimación a priori, no resultado).</li>
    <li><b>Especulativa:</b> con ese potencial coulombiano el amortiguamiento cerca del umbral tendería a un valor finito (factor de Sommerfeld) en lugar de anularse como en la sección 3.</li>
  </ol>
</section>'''

page = f'''<title>Umbral axiónico radial</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
/* Layout: una columna de lectura (~68ch) con figuras algo más anchas; misma identidad que el informe principal */
:root {{
  --paper: #f4f6f9; --sheet: #ffffff; --ink: #161e2e; --muted: #586378; --rule: #d8dde6;
  --accent: #2b4bb3; --ok: #1d7451; --ok-bg: #e4f3eb; --flag: #9a5400; --flag-bg: #fdf0dc;
  --bad: #a3262f; --bad-bg: #fbe5e6; --code-bg: #eaeef4;
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --grid: #e6e9ef;
  --f-display: "Newsreader", "Iowan Old Style", Georgia, serif;
  --f-body: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  --f-mono: "IBM Plex Mono", ui-monospace, "SFMono-Regular", Menlo, monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --paper: #0d121b; --sheet: #141b27; --ink: #e3e8f1; --muted: #98a3b7; --rule: #263042;
    --accent: #93a8ff; --ok: #62c99c; --ok-bg: #12291f; --flag: #f1a84d; --flag-bg: #2b2112;
    --bad: #f08a90; --bad-bg: #2d1518; --code-bg: #1b2433;
    --s1: #3987e5; --s2: #d95926; --s3: #199e70; --grid: #212a3a; color-scheme: dark;
  }}
}}
:root[data-theme="dark"] {{
  --paper: #0d121b; --sheet: #141b27; --ink: #e3e8f1; --muted: #98a3b7; --rule: #263042;
  --accent: #93a8ff; --ok: #62c99c; --ok-bg: #12291f; --flag: #f1a84d; --flag-bg: #2b2112;
  --bad: #f08a90; --bad-bg: #2d1518; --code-bg: #1b2433;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --grid: #212a3a; color-scheme: dark;
}}
body {{ background: var(--paper); color: var(--ink); font: 16px/1.6 var(--f-body); }}
.wrap {{ max-width: 860px; margin: 0 auto; padding-left: max(16px, env(safe-area-inset-left, 0px));
  padding-right: max(16px, env(safe-area-inset-right, 0px)); padding-block: 40px 64px; display: grid; gap: 28px; }}
header {{ display: grid; gap: 10px; }}
.eyebrow {{ font: 500 12px/1 var(--f-mono); letter-spacing: .08em; text-transform: uppercase; color: var(--muted); }}
h1 {{ font: 600 clamp(30px, 5vw, 44px)/1.1 var(--f-display); margin: 0; text-wrap: balance; }}
h2 {{ font: 600 25px/1.2 var(--f-display); margin: 0; text-wrap: balance; }}
h3 {{ font: 600 16px/1.3 var(--f-body); margin: 0; }}
.lede {{ font-size: 18px; color: var(--muted); max-width: 68ch; margin: 0; }}
section {{ display: grid; gap: 14px; }}
section > p, section > ul, section > ol {{ max-width: 68ch; margin: 0; }}
.sheet {{ background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; padding: 20px 22px; display: grid; gap: 10px; min-width: 0; }}
.summary ul {{ margin: 0; padding-left: 20px; display: grid; gap: 6px; }}
code, .mono {{ font-family: var(--f-mono); font-size: .92em; }}
code {{ background: var(--code-bg); padding: 1px 5px; border-radius: 3px; }}
figure {{ margin: 0; display: grid; gap: 8px; min-width: 0; }}
figure svg {{ width: 100%; height: auto; display: block; background: var(--sheet); border: 1px solid var(--rule); border-radius: 6px; }}
figcaption {{ font-size: 14px; color: var(--muted); max-width: 72ch; }}
.legend {{ display: flex; flex-wrap: wrap; gap: 6px 18px; font-size: 13px; color: var(--muted); }}
.legend span {{ display: inline-flex; align-items: center; gap: 6px; }}
.sw {{ width: 10px; height: 10px; border-radius: 50%; display: inline-block; }}
.sw.s1 {{ background: var(--s1); }} .sw.s2 {{ background: var(--s2); }} .sw.s3 {{ background: var(--s3); }}
.sw.ref {{ background: transparent; border: 1.5px solid var(--muted); }}
svg .grid {{ stroke: var(--grid); stroke-width: 1; }}
svg .frame {{ fill: none; stroke: var(--rule); }}
svg .zero {{ stroke: var(--muted); stroke-width: 1; opacity: .6; }}
svg .thr {{ stroke: var(--accent); stroke-width: 1.5; stroke-dasharray: 5 4; }}
svg .unst {{ stroke: var(--bad); stroke-width: 1.5; stroke-dasharray: 2 3; }}
svg .tick {{ fill: var(--muted); font: 12px var(--f-mono); }}
svg .axlabel {{ fill: var(--ink); font: 500 13px var(--f-body); }}
svg .note {{ fill: var(--ink); font: 12.5px var(--f-body); }}
svg .zone {{ fill: var(--muted); font: italic 12px var(--f-body); }}
svg .s1-line {{ fill: none; stroke: var(--s1); stroke-width: 2; }}
svg .ref-line {{ fill: none; stroke: var(--muted); stroke-width: 1.5; stroke-dasharray: 6 4; }}
svg .s1-dot {{ fill: var(--s1); stroke: var(--sheet); stroke-width: 1.5; }}
svg .s2-dot {{ fill: var(--s2); stroke: var(--sheet); stroke-width: 1.5; }}
svg .s3-dot {{ fill: var(--s3); stroke: var(--sheet); stroke-width: 1.5; }}
svg .s1-ring {{ fill: var(--sheet); stroke: var(--s1); stroke-width: 2.5; }}
svg .s2-ring {{ fill: var(--sheet); stroke: var(--s2); stroke-width: 2.5; }}
svg .ref-dot {{ fill: var(--sheet); stroke: var(--muted); stroke-width: 1.5; }}
svg circle:hover {{ stroke: var(--ink); stroke-width: 2; }}
.tbl {{ overflow-x: auto; }}
table {{ border-collapse: collapse; font-size: 14px; font-variant-numeric: tabular-nums; min-width: 420px; }}
th, td {{ text-align: left; padding: 6px 14px 6px 0; border-bottom: 1px solid var(--rule); }}
th {{ font-weight: 600; color: var(--muted); font-size: 12px; letter-spacing: .04em; text-transform: uppercase; }}
.callout {{ border-left: 3px solid var(--flag); background: var(--flag-bg); padding: 14px 18px; border-radius: 0 6px 6px 0; display: grid; gap: 8px; }}
.callout.neutral {{ border-left-color: var(--muted); background: var(--code-bg); }}
.callout p {{ margin: 0; max-width: 68ch; }}
.checks {{ list-style: none; margin: 0; padding: 0; display: grid; gap: 8px; }}
.checks li {{ display: grid; grid-template-columns: auto 1fr; gap: 10px; align-items: start; font-size: 14px; }}
.checks li > span:last-child {{ min-width: 0; overflow-wrap: anywhere; }}
.pill {{ font: 600 11px/1 var(--f-mono); padding: 4px 6px; border-radius: 3px; letter-spacing: .05em; }}
.ok .pill {{ background: var(--ok-bg); color: var(--ok); }}
.bad .pill {{ background: var(--bad-bg); color: var(--bad); }}
.status {{ font: 600 12px var(--f-mono); color: var(--ok); }}
.status.pend {{ color: var(--flag); }}
footer {{ font-size: 13px; color: var(--muted); border-top: 1px solid var(--rule); padding-top: 14px; }}
a {{ color: var(--accent); }}
:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
</style>
<div class="wrap">
<header>
  <div class="eyebrow">Workshop IA × física · pasos 14–15 · modelos radiales del Paper I</div>
  <h1>Del modo cuasi-normal al modo ligado: el umbral axiónico en las oscilaciones radiales</h1>
  <p class="lede">Respuesta a la pregunta de José A. Pons: cómo pasa un modo radial de cuasi-normal a evanescente cuando su frecuencia cruza la masa del axión, y cómo se interpretan las distintas ramas del espectro. Hipótesis anotadas antes de calcular; todas las cifras salen de <code>results/radial_branches.json</code>.</p>
</header>

<section class="sheet summary">
  <h3>Resumen</h3>
  <ul>
    <li>Para <code>l = 0</code> no hay radiación gravitatoria, así que el único canal de pérdida es la emisión de axiones, con <code>k_out² = ω² − m_a²</code>. Escrita en la variable <code>k = k_out</code>, la función de dispersión del modelo (Ec. 4.13) es univaluada y no hay que elegir hoja.</li>
    <li>La transición sigue la secuencia típica de onda con umbral: el par cuasi-normal <code>k = ±a − ib</code> colisiona en el eje imaginario negativo (entre <code>x₀ = {fnum(x_prev, 4)}</code> y <code>{fnum(x_coll, 4)}</code>) y se separa en un <b>estado ligado</b> (<code>Im k &gt; 0</code>, frecuencia real) y un <b>estado virtual</b> (<code>Im k &lt; 0</code>, frecuencia real pero no físico).</li>
    <li>Cerca del umbral, <code>|Im w| ∝ (Re w − 1)^{fnum(s_fit, 2)}</code>, próximo a ½. El amortiguamiento no se anula exactamente en <code>ω = m_a</code>, sino a una distancia <code>∼ (Im k)²</code>.</li>
    <li>Las ramas fluida y axiónica no se cruzan: la separación mínima es lineal en el acoplo, con <code>Δw²_min/λ = {fnum(min(ratios), 2)}–{fnum(max(ratios), 2)}</code> para λ² = 10⁻⁶…10⁻³.</li>
    <li>Hallazgo no previsto: linealizado en <code>θ = 0</code>, el modelo contiene un modo inestable (<code>w² = {fnum(w2u, 2)}</code> con los parámetros de la figura), que no aparece en la ventana <code>Im w ≤ 0</code> habitual. Lo planteamos como pregunta abajo.</li>{star_summary}
  </ul>
</section>

<section>
  <h2>1. Planteamiento</h2>
  <p>Usamos el modelo simplificado de la Sec. IV del Paper I (Ec. 4.13, ya reproducido en el paso 4): fondo homogéneo, simetría plana, sin perturbación métrica, con <code>w = ω/m_a</code> y <code>x₀ = m_a R</code>. Los parámetros son los de la Fig. 3: <code>α = {fnum(alpha, 0)}</code>, <code>c_s² = {fnum(D["params"]["cs2"], 1)}</code> y <code>λ² = 10⁻³</code>. Para mover un modo fluido a través del umbral variamos <code>x₀</code>, lo que equivale a cambiar <code>m_a R</code>.</p>
  <p>La función depende de <code>ω</code> solo a través de <code>w²</code>, salvo el término de onda exterior <code>i k_out</code>. Con <code>w² = 1 + k²</code> queda como función entera de <code>k</code>, y cada raíz se clasifica sin ambigüedad:</p>
  <div class="tbl"><table>
    <tr><th>Posición de la raíz en el plano k</th><th>Onda exterior</th><th>Interpretación</th></tr>
    <tr><td><code>Re k &gt; 0</code>, <code>Im k &lt; 0</code></td><td>saliente, crece con r</td><td>modo cuasi-normal (amortiguado)</td></tr>
    <tr><td><code>Re k = 0</code>, <code>Im k &gt; 0</code></td><td>evanescente <code>e^(−κr)</code></td><td>modo ligado: frecuencia real, normalizable</td></tr>
    <tr><td><code>Re k = 0</code>, <code>Im k &lt; 0</code></td><td>crece <code>e^(+κr)</code></td><td>estado virtual: frecuencia real, no es un modo físico</td></tr>
  </table></div>
  <p>En la variable <code>w</code> esa distinción se pierde: un estado ligado y uno virtual tienen ambos <code>Im w = 0</code>. La expresión <code>w√(1 − 1/w²)</code> del paso 4 pone el corte en el segmento real <code>(−1, 1)</code>. Por debajo del umbral y en el semiplano inferior da la hoja no física. Nuestro primer seguimiento de raíces cayó precisamente en el estado virtual.</p>
</section>

<section>
  <h2>2. La transición</h2>
  <figure>
    <div class="legend"><span><i class="sw s1"></i>par cuasi-normal ±a − ib</span><span><i class="sw s3"></i>estado ligado</span><span><i class="sw s2"></i>estado virtual</span></div>
    {fig1}
    <figcaption>Raíces del modo fluido fundamental en el plano k, para x₀ de {fnum(pair[0][0], 2)} a {fnum(after[-1]["x0"], 2)} (λ² = 10⁻³). El par baja hacia el eje con <code>b</code> casi constante ({fnum(b_lo, 4)}–{fnum(b_hi, 4)}) mientras <code>Re k → 0</code>, y colisiona entre <code>x₀ = {fnum(x_prev, 4)}</code> y <code>{fnum(x_coll, 4)}</code>. Después una raíz sube por <code>k = 0</code> (ligado) y la otra baja (virtual). Pasa el cursor por los puntos para ver los valores.</figcaption>
  </figure>
  <p>En <code>x₀ = 0.86</code> el estado ligado está en <code>κ = +{fnum(bound86, 4)}</code> (<code>w = {fnum(w86, 5)}</code>, real a precisión de máquina) y el virtual en <code>κ = {fnum(virt86, 4)}</code>. Sin acoplo (λ = 0), cada modo fluido aparece como un par simétrico <code>±κ</code>, porque el fluido no ve el exterior. Es el acoplo con el axión lo que rompe esa simetría y fija qué raíz es la física.</p>
</section>

<section>
  <h2>3. El amortiguamiento cerca del umbral</h2>
  <figure>
    {fig3}
    <figcaption>Tasa de amortiguamiento frente a la distancia al umbral (λ² = {D["params"]["lam2_damping"]:.0e}), en la zona asintótica <code>Re k &gt; 5|Im k|</code>. Ajuste: exponente {fnum(s_fit, 3)}; la línea discontinua es la pendiente ½ como referencia.</figcaption>
  </figure>
  <p>La identidad exacta <code>Im w = −Re k · Im k / Re w</code> explica el comportamiento: con <code>Im k</code> casi constante, la tasa es proporcional a <code>k_out ∝ √(ω² − m_a²)</code>, el espacio de fases de la onda emitida. El exponente ajustado queda algo por debajo de ½ porque <code>Im k</code> varía un ~18 % en el tramo. La ley se rompe cuando <code>Re k ∼ |Im k|</code>: ahí se produce la colisión, a una distancia del umbral de orden <code>(Im k)²</code>, minúscula con acoplo débil.</p>
</section>

<section>
  <h2>4. Ramas del espectro y cruces evitados</h2>
  <figure>
    <div class="legend"><span><i class="sw ref"></i>sin acoplo (λ = 0)</span><span><i class="sw s1"></i>con acoplo (λ² = 10⁻³)</span></div>
    {fig2}
    <figcaption>Estados ligados (<code>Im k &gt; 0</code>) frente a <code>x₀</code>. Por debajo de la línea roja la frecuencia es imaginaria (modo que crece). No se muestra la rama fundamental del pozo axiónico, con w² entre {fnum(max(deep), 1)} y {fnum(min(deep), 1)} en todo el intervalo (sección 5).</figcaption>
  </figure>
  <p>Hay dos familias de modos:</p>
  <ul>
    <li>La <b>rama fluida</b> sigue <code>w ≈ c_s(n + ½)π/x₀</code>, con amortiguamiento de orden λ² por encima del umbral.</li>
    <li>La <b>rama axiónica</b> es la escalera de estados del pozo interior <code>m²_in = 1 − α &lt; 0</code>. El segundo estado se liga en el umbral cuando <code>√α · x₀ = π</code> (<code>x₀ = {fnum(math.pi / math.sqrt(alpha), 4)}</code>), baja y cruza <code>w² = 0</code> poco después.</li>
  </ul>
  <p>Sin acoplo, las ramas se cruzan. Con acoplo se repelen e intercambian carácter, y la separación mínima escala linealmente con λ:</p>
  <div class="tbl"><table>
    <tr><th>λ²</th><th>Δw²<sub>min</sub></th><th>Δw²<sub>min</sub>/λ</th><th>x₀ del cruce</th></tr>
    {gap_rows}
  </table></div>
</section>

<section>
  <h2>5. Un modo inestable en el modelo linealizado en θ = 0</h2>
  <p>Con los parámetros de la Fig. 3 (<code>x₀ = 0.6</code>) hay un estado ligado con <code>w² = {fnum(w2u, 3)}</code>, es decir, una perturbación que crece como <code>e^({fnum(grow, 2)} m_a t)</code>. Es el estado fundamental del pozo taquiónico interior, con onda oscilante dentro y evanescente fuera. Aparece para <code>x₀ &gt; {fnum(xc, 5)}</code>, valor que coincide con la condición analítica <code>arctan(1/√(α−1))/√(α−1) = {fnum(xc_an, 5)}</code>.</p>
  <p>Nuestro paso 4 no lo detectaba porque buscaba raíces solo en <code>Im w ≤ 0.01</code>. Linealizando en <code>θ = π</code> (masa efectiva interior positiva, sin acoplo lineal) desaparece: no hay estados ligados ni inestables.</p>
  <div class="callout">
    <h3>Pregunta para el grupo</h3>
    <p>¿El modelo simplificado se pensó como perturbación del estado sin condensar? En ese caso, el modo inestable es la inestabilidad que produce la condensación y forma parte de la historia. Si se pensó como perturbación del condensado, ¿no debería linealizarse en <code>θ = π</code>? Allí <code>∂m*/∂a = 0</code> y <code>m²_eff = m_a²(α − 1)/√(1 − β) &gt; 0</code> (ya lo apuntamos al reproducir el paso 3).</p>
  </div>
</section>

{star_html if star_html else PENDING}

<section>
  <h2>7. Verificación del modelo simplificado</h2>
  <p><b>{n_ok}/{len(checks)}</b> comprobaciones superadas en <code>axion_ns/step14_radial_branches.py</code> (unos 30 s, incluido en <code>run_all</code>). Los dos controles negativos cambian la física (onda entrante; signo de la masa interior) y deben romper el resultado.</p>
  <ul class="checks">{check_rows}</ul>
  <div class="callout neutral">
    <p><b>Límites.</b> La función de dispersión es la Ec. (4.13) tal como la transcribimos en el paso 4, no una compilación desde SymPy. Las comprobaciones de la estrella realista están al final de la sección 6.</p>
  </div>
</section>

<section>
  <div class="callout neutral">
    <h3>Nota sobre el informe anterior de «tres líneas de investigación»</h3>
    <p>Lo hemos retirado. Dos revisiones independientes y una auditoría del código mostraron que los modos g, la rotación y el puente axión-magnetar se apoyaban en constantes ajustadas a mano y en premisas físicas incorrectas: un salto de <code>m_N*</code> diez veces mayor, <code>f_GW = 2ν</code> en la resonancia de marea y <code>μ₅</code> a partir de <code>∇a</code> en vez de <code>∂ₜa</code>. Las comprobaciones eran circulares. El análisis completo está en el cuaderno de laboratorio, y desde entonces hemos endurecido las reglas de verificación. Los pasos 1–10 (informe principal) no se ven afectados.</p>
  </div>
</section>

<footer>Generado por <code>report/build_radial.py</code> desde <code>results/radial_branches.json</code> (paso 14) y <code>results/radial_star.json</code> (paso 15). Referencia: A. Gómez-Bañón, P. Pnigouras y J. A. Pons, <a href="https://arxiv.org/abs/2510.16245">arXiv:2510.16245</a>.</footer>
</div>
'''

for out in (ROOT / 'report' / 'informe_radial.html', ROOT / 'results' / 'informe_radial.html'):
    out.write_text(page)
print('informe_radial.html:', len(page) // 1024, 'KB;', n_ok, '/', len(checks), 'checks')
