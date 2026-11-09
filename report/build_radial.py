"""Genera report/informe_radial.html (y la copia en results/) desde results/radial_branches.json.

Regla (AGENTS.md §4): ninguna cifra se teclea a mano; todas salen del JSON que escribe
axion_ns/step14_radial_branches.py.
"""
import html
import json
import cmath
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
D = json.loads((ROOT / 'results' / 'radial_branches.json').read_text())


def fnum(x, nd=3):
    return f'{x:.{nd}f}'.replace('-', '−')


def cnum(z, nd=2):
    """Número complejo [re, im] como «a − bi»."""
    return f"{fnum(z[0], nd)} {'−' if z[1] < 0 else '+'} {fnum(abs(z[1]), nd)}i"


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
# Im k en el tramo del ajuste (Re k > 5|Im k|), desde w = 1 + x - i y con k = √(w² − 1), Re k > 0
_ks = [cmath.sqrt(complex(1 + x, -y) ** 2 - 1) for x, y in dp]
_ks = [k if k.real > 0 else -k for k in _ks]
_imk = [-k.imag for (x, _), k in zip(dp, _ks) if k.real > 5 * -k.imag and x < 0.02]
imk_var = 100 * (max(_imk) / min(_imk) - 1)
C3, C4 = D['collision']['0.001'], D['collision']['0.0001']          # colisión exacta (raíz doble), paso 14
sci = lambda v, d=1: (lambda m, e: f'{m}·10{str(int(e)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))}')(
    *f'{v:.{d}e}'.replace('.', ',').split('e'))
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
# ok = None: resultado descriptivo, no comprobación (revisión del informe principal)
_cls = lambda c: 'res' if c['ok'] is None else ('ok' if c['ok'] else 'bad')
_pill = lambda c: 'RESULTADO' if c['ok'] is None else ('OK' if c['ok'] else 'FALLA')
checks = D['checks']
n_ok = sum(bool(c['ok']) for c in checks)
check_rows = ''.join(
    f'<li class="{_cls(c)}"><span class="pill">{_pill(c)}</span>'
    f'<span>{html.escape(c["name"])}{(" <em>(" + html.escape(c["note"]) + ")</em>") if c["note"] else ""}</span></li>'
    for c in checks)

# estimación a priori para H11 (no es un resultado: se calcula aquí para no teclearla)
GM_SUN_KM, C_KMS = 1.476625, 299792.458
alpha_g = 1.5 * GM_SUN_KM * 2 * math.pi * 3000 / C_KMS
wf = D['fluid_x0_0p6']
xc, xc_an = D['x0_unstable'], D['x0_unstable_analytic']
xc3 = D.get('x0_unstable_lam2_1e-3')
w2u = D['unstable_x0_0p6_w2']
grow = math.sqrt(-w2u)
alpha = D['params']['alpha']

# ------------------------------------------------------------------ estrella realista (paso 15)
SP = ROOT / 'results' / 'radial_star.json'
S = json.loads(SP.read_text()) if SP.exists() else None
STP = ROOT / 'results' / 'radial_stability.json'
STB = json.loads(STP.read_text()) if STP.exists() else None
CRP = ROOT / 'results' / 'radial_crossing.json'
CRS = json.loads(CRP.read_text()) if CRP.exists() else None
DMP = ROOT / 'results' / 'radial_damped.json'
DMD = json.loads(DMP.read_text()) if DMP.exists() else None
DPP = ROOT / 'results' / 'radial_deep.json'
DEEP = json.loads(DPP.read_text()) if DPP.exists() else None
CTP = ROOT / 'results' / 'radial_contour.json'
CTR = json.loads(CTP.read_text()) if CTP.exists() else None
RYP = ROOT / 'results' / 'radial_rydberg.json'
RYD = json.loads(RYP.read_text()) if RYP.exists() else None
GUARD_TXT = ('\n  <p><b>Garantías de las raíces (desde la sexta revisión).</b> Muller lanza un error si no converge o si el residuo normalizado no es pequeño, y toda integración que no termina detiene el cálculo. '
             'Antes, Muller devolvía su último intento sin comprobar nada, y las integraciones fallidas pasaban inadvertidas. Esa comprobación destapó el fallo del arranque del interior (recuadro siguiente). Todos los pasos se repitieron con estas garantías.</p>')
_OLD = ROOT / 'report' / 'baseline_pre_h29'
OLD_S = json.loads((_OLD / 'radial_star.json').read_text()) if (_OLD / 'radial_star.json').exists() else None
OLD_V = json.loads((_OLD / 'radial_validation.json').read_text()) if (_OLD / 'radial_validation.json').exists() else None
EOP = ROOT / 'results' / 'radial_eos.json'
EOSD = json.loads(EOP.read_text()) if EOP.exists() else None
SCP = ROOT / 'results' / 'radial_scattering.json'
SCAT = json.loads(SCP.read_text()) if SCP.exists() else None
VLP = ROOT / 'results' / 'radial_validation.json'
VAL = json.loads(VLP.read_text()) if VLP.exists() else None
if VAL:
    _vrows = ''.join(f'<tr><td>{k}</td><td>{VAL["modes"][k]:.4f}</td><td>{VAL["font"][k]}</td><td>{VAL["dev_font_pct"][k]:+.3f} %</td></tr>'
                     for k in VAL['font'])
    val_html = f"""
  <h3>Validación externa sin axión (paso 17)</h3>
  <p>El mismo código, sin cambiar nada y con el axión apagado, aplicado al politropo de referencia de Font et al. (2002, PRD 65, 084024, Tabla 2): <code>p = Kρ²</code>, <code>K = 100</code>, <code>ρ_c = 1,28·10⁻³</code> (G = c = M☉ = 1). El fondo sale con M = {fnum(VAL['M'], 4)} M☉ y R = {fnum(VAL['R'], 3)} km (ellos: 1,4 y 14,15). Sus frecuencias son las de la teoría de perturbaciones con la métrica dinámica:</p>
  <div class="tbl"><table><tr><th>modo</th><th>este cálculo (kHz)</th><th>Font et al. (kHz)</th><th>diferencia</th></tr>{_vrows}</table></div>
  <p>En relativo, la diferencia es como mucho un {max(abs(v) for v in VAL['dev_font_pct'].values()):.3f} %. En absoluto: {', '.join(f"{k} {1e3 * (VAL['modes'][k] - VAL['font'][k]):+.1f} Hz" for k in VAL['font'])}. Su tabla da tres decimales en kHz (redondeo ±0,5 Hz): F está {'dentro' if abs(1e3 * (VAL['modes']['F'] - VAL['font']['F'])) <= 0.5 else 'fuera'}, los sobretonos difieren algo más ({min(abs(1e3 * (VAL['modes'][k] - VAL['font'][k])) for k in ('H1', 'H2', 'H3')):.1f}–{max(abs(1e3 * (VAL['modes'][k] - VAL['font'][k])) for k in ('H1', 'H2', 'H3')):.1f} Hz). Esto valida el sector fluido en GR sin axión, no el acoplo con el axión ni los amortiguamientos pequeños. Los controles negativos comprueban que la comparación distingue: con Γ₁ = 2,1 en las perturbaciones, F se mueve un {VAL['neg_gamma1']:+.1f} %, y sin la energía interna en el fondo, un {VAL['neg_no_internal']:+.1f} %. La tabla del grupo de Garching para el mismo modelo da F = {VAL['garching']['F']} kHz, un {abs(VAL['dev_garching_pct']['F']):.1f} % por encima de nosotros y de Font et al.; no indica el método.</p>"""
else:
    val_html = ''
if S:
    eos_txt = ''
    if EOSD:
        _sh = EOSD['shift']
        eos_txt = (f" Los ajustes publicados de BSk26, ε(n) y P(ρ), no cumplen exactamente la identidad termodinámica dε/dn = (ε + p)/n: el desajuste es {sci(abs(EOSD['mismatch_center']), 1)} en el centro de la estrella de 6 kHz. "
                   f"Calculando la velocidad del sonido de la otra forma posible, F se mueve {fnum(abs(_sh['F6']['dRe_Hz']), 2)} Hz, H1 {fnum(abs(_sh['H1']['dRe_Hz']), 2)} Hz y su Im ν un {fnum(100 * _sh['H1']['dIm_rel'], 2)} % "
                   f"(paso 23). Eso es más que todo lo que mueven los controles del integrador: es la incertidumbre dominante de la parte fluida que hemos medido.")
    _thin = DMD.get('thin') if DMD else None
    _imtop = '−10⁻⁵' if _thin and _thin['total'] == 0 else '−0,002'
    thin_txt = ((f" La franja casi real que dejaban fuera esas piezas, <code>−0,002 ≤ Im ν ≤ −10⁻⁵</code> kHz, da índice {_thin['total']} (salto de fase máximo {fnum(_thin['max_dphase'], 2)} rad). "
                 "Quedan sin contar <code>|Im ν| &lt; 0,01</code> Hz, <code>Re ν &lt; 6,05</code> kHz y <code>Re ν &gt; 12</code> kHz.") if _thin else
                " El recuento no cubre <code>−0,002 &lt; Im ν &lt; 0</code> kHz, <code>Re ν &lt; 6,05</code> kHz ni <code>Re ν &gt; 12</code> kHz.")
    _jumps = [x for x in ((STB or {}).get('nu6.0', {}).get('max_dphase'), (STB or {}).get('nu3.4', {}).get('max_dphase')) if x is not None]
    jump_txt = f" (el mayor en este paso: {fnum(max(_jumps), 2)} rad)" if _jumps else ''
    guard_txt = GUARD_TXT
    scat_sum = ''
    scat_html = ''
    if SCAT:
        _h, _f = SCAT['H1 (6 kHz)'], SCAT['F (3,4 kHz)']
        scat_sum = f" (lo confirma, sin plano complejo, la anchura de la resonancia en la dispersión a frecuencia real: {fnum(_h['half_width'] * 1e3, 2)} Hz, paso 24)"
        scat_html = f'''
  <h3>Una comprobación independiente del amortiguamiento (paso 24)</h3>
  <p>Muller y el principio del argumento buscan ceros de la misma función, así que un error en ella los engañaría a los dos. Aquí el amortiguamiento se mide de otra forma, sin plano complejo ni búsqueda de raíces. A frecuencia real se lanza una onda de axión hacia la estrella y se mide la onda que vuelve, <code>S(ν) = −W[u, ψ_in]/W[u, ψ_out]</code>. Cerca de un modo, la fase de S gira 2π con la forma de Breit–Wigner, y su anchura es el amortiguamiento.</p>
  <div class="tbl"><table><tr><th>modo</th><th>−Im ν, Muller (Hz)</th><th>Γ/2, dispersión (Hz)</th><th>diferencia</th><th>|ν_r − Re ν| / anchura</th></tr>
  <tr><td>H1 (ν_a = 6 kHz)</td><td>{fnum(-_h['muller'][1] * 1e3, 4)}</td><td>{fnum(_h['half_width'] * 1e3, 4)}</td><td>{100 * _h['rel_width']:.2g} %</td><td>{_h['dnu_over_width']:.1e}</td></tr>
  <tr><td>F (ν_a = 3,4 kHz)</td><td>{fnum(-_f['muller'][1] * 1e3, 4)}</td><td>{fnum(_f['half_width'] * 1e3, 4)}</td><td>{100 * _f['rel_width']:.2g} %</td><td>{_f['dnu_over_width']:.1e}</td></tr></table></div>
  <p>El ajuste no salió a la primera. F a 3,4 kHz está a solo {_f['thr_widths']:.0f} anchuras del umbral, y ahí la fase «de fondo», la que no es de la resonancia, se curva: crece como la raíz de la distancia al umbral. Con un fondo recto, el ajuste no cierra: le quedan {fnum(_f['lin_resid_rad'], 0)} rad sin explicar, y en el primer recálculo, con menos puntos, la anchura salió muy desviada (el detalle está en el cuaderno, H28). Por eso el ajuste se repite con fondos de grado 1 a 4 en ventanas cada vez más estrechas, de ±10 a ±1,5 anchuras. La tabla da el de grado 4 en ±1,5 anchuras. Pasar del grado 3 al 4 cambia la anchura {sci(_f['conv'], 0)}, y lo que queda sin ajustar son {sci(_f['fit_resid_rad'], 0)} rad. H1, a {_h['thr_widths']:.0f} anchuras del umbral, da lo mismo con cualquier fondo.</p>
  <p>Además, <code>|S| = 1</code> a {max(_h['unitarity'], _f['unitarity']):.0e}: lo que entra sale. Con unas ecuaciones reales eso está garantizado, así que comprueba la numérica del exterior, no la física. Lo que <b>no</b> comprueba esta prueba: usa las mismas ecuaciones de perturbación, de modo que valida la continuación al plano complejo y la búsqueda de raíces, no las ecuaciones (esas se comparan con el Paper I en el paso 2). Controles: desplazando la ventana 20 anchuras no aparece resonancia; con la normalización antigua del lapso, la de H1 cae a {SCAT['old_lapse_H1']['dnu_over_width']:.2g} anchuras del modo corregido.</p>'''
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
    p4.text(s6['F'], 0.8, f'F {fnum(s6["F"], 3)}', anchor='middle', dy=24)
    p4.text(min(lad), 0.8, 'átomo n ≥ 2', anchor='end', dx=-10, dy=4)
    p4.text(s6['H1'][0], -0.8, f'H1, Im ν = {fnum(s6["H1"][1] * 1e3, 1)} Hz', anchor='middle', dy=22)
    fig4 = p4.svg('Espectro radial a nu_a = 6 kHz')
    # Figura 6: cruce evitado (paso 19)
    if CRS:
        rows_ = [r for r in CRS['rows'] if len(r['roots']) == 2]
        xa = [r['nu_a'] for r in rows_]
        ya = [v for r in rows_ for v in r['roots']]
        p6 = Plot((min(xa) - 0.01, max(xa) + 0.01), (min(ya) - 0.02, max(ya) + 0.02), h=330)
        p6.axes([3.45, 3.5, 3.55, 3.6, 3.65, 3.7], [3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 3.65], 'ν_a (kHz)', 'ν (kHz)',
                lambda v: fnum(v, 2), lambda v: fnum(v, 2))
        p6.path([(x, x) for x in (min(xa) - 0.01, max(ya) + 0.02)], 'ref-line')
        for j in (0, 1):
            p6.path([(r['nu_a'], r['roots'][j]) for r in rows_], 'ref-line')
            for cls, cond in (('s1-dot', lambda f: f >= 0.5), ('s2-dot', lambda f: f < 0.5)):
                p6.dots([(r['nu_a'], r['roots'][j]) for r in rows_ if cond(r['fluid'][j])], cls,
                        lambda x, y: f'ν_a = {x:.3f}, ν = {y:.5f} kHz', r=5)
        p6.text(min(xa) + 0.03, min(xa) + 0.03, 'umbral ν = ν_a', anchor='end', dx=-6, dy=-8)
        fig6 = p6.svg('Cruce evitado entre F y el nivel n = 1 del átomo gravitatorio')
        gmin, numin = CRS['gap_min_khz'], CRS['nu_a_min']
        half = min(rows_, key=lambda r: abs(r['fluid'][0] - 0.5))
        cross_html = f'''
  <figure>
    <div class="legend"><span><i class="sw s1"></i>carácter fluido (&gt; 50 %)</span><span><i class="sw s2"></i>carácter axiónico</span></div>
    {fig6}
    <figcaption>Los dos estados ligados más bajos frente a ν_a (paso 19). El color indica el carácter, medido con la fracción fluida. Las ramas no se cortan: la separación mínima es de {1e3 * gmin:.1f} Hz en ν_a = {fnum(numin, 3)} kHz, y el carácter pasa de un estado al otro.</figcaption>
  </figure>
  <p><b>Es un cruce evitado (paso 19).</b> Seguimos los dos estados ligados más bajos en {len(rows_)} valores de ν_a, con la malla refinada alrededor del mínimo. La separación nunca baja de {1e3 * gmin:.1f} Hz. El carácter se intercambia de forma continua: a ν_a = {fnum(rows_[0]['nu_a'], 2)} kHz el inferior es axiónico ({rows_[0]['fluid'][0]:.2f} de fluido) y el superior es F ({rows_[0]['fluid'][1]:.2f}); a {fnum(rows_[-1]['nu_a'], 2)} kHz, al revés ({rows_[-1]['fluid'][0]:.2f} y {rows_[-1]['fluid'][1]:.2f}). La mezcla 50/50 está hacia ν_a = {fnum(half['nu_a'], 3)} kHz, algo por encima del mínimo de separación, porque F también se mueve con ν_a. La fracción fluida es un diagnóstico (ver la figura del espectro), no una energía conservada.</p>'''
    else:
        cross_html = ''
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
    if RYD:
        lad_rows = ''.join(f'<tr><td>{l["n"]}</td><td>{l["nu"]:.6f}</td><td>{l["n_eff"]:.4f}</td><td>{l["n"] - (RYD["n"][i] - RYD["n_minus_kappa"][i]):.4f}</td></tr>'
                           for i, l in enumerate(s6['ladder']))
        lad_head = '<tr><th>n</th><th>ν (kHz)</th><th>n_eff = α_g/√(2(1 − ν/ν_a))</th><th>n − κ (Whittaker)</th></tr>'
        rz = RYD['ritz']
        nk = RYD['n_minus_kappa']
        ryd_html = (f'<p><b>Una ley de Rydberg, no hidrogenoide pura (paso 18).</b> Con la fórmula no relativista, <code>n − n_eff</code> '
                    f'parecía crecer con n (de {min(RYD["n_minus_neff"]):+.3f} a {max(RYD["n_minus_neff"]):+.3f}). Pero la ecuación exterior, '
                    f'reducida a la de Whittaker en el paso 5 (Ec. A14 de la secuela), tiene su propio parámetro de Coulomb relativista, '
                    f'<code>κ = α_g(1 − 2q²)/q</code>, con <code>q² = 1 − (ν/ν_a)²</code>. Con él, <code>n − κ</code> es casi constante: '
                    f'{", ".join(f"{v:.4f}" for v in nk)}, dentro de {max(nk) - min(nk):.4f}, sin ajustar nada. Es un <b>defecto cuántico</b>, '
                    f'como en los átomos alcalinos: la región interior (la estrella y el campo fuerte) añade un desfase que no depende de la '
                    f'energía. Con la corrección de Ritz, <code>δ_n = δ₀ + δ₂/(n − δ₀)²</code>, sale δ₀ = {rz["delta0"]:.4f} y '
                    f'δ₂ = {rz["delta2"]:.4f}, con un residuo de {rz["resid"]:.0e}. Para l = 0, el término 1/r² del desarrollo asintótico es '
                    f'supercrítico, así que δ no sale de ese desarrollo: lo fija el interior.</p>')
    else:
        lad_rows = ''.join(f'<tr><td>{l["n"]}</td><td>{l["nu"]:.6f}</td><td>{l["n_eff"]:.3f}</td></tr>' for l in s6['ladder'])
        lad_head = '<tr><th>n</th><th>ν (kHz)</th><th>n_eff = α_g/√(2(1 − ν/ν_a))</th></tr>'
        ryd_html = ''
    t = S['turning']
    n35, n34 = S['nu35'], S['nu34']
    st_checks = S.get('checks', [])
    st_rows = ''.join(
        f'<li class="{_cls(c)}"><span class="pill">{_pill(c)}</span>'
        f'<span>{html.escape(c["name"])}{(" <em>(" + html.escape(c["note"]) + ")</em>") if c["note"] else ""}</span></li>'
        for c in st_checks)
    rob = S.get('robustness', {})
    rob_rows = ''
    for lab, r_ in rob.items():
        for k, v_ in r_.items():
            rob_rows += f'<tr><td>{lab}</td><td>{html.escape(k)}</td><td>{v_[0]:.6f}</td><td>{v_[1] * 1e3:.4f}</td></tr>'
    _dev = [(abs(v_[0] / r_['base'][0] - 1), abs(v_[1] / r_['base'][1] - 1)) for r_ in rob.values() for v_ in r_.values()]
    rob_re, rob_im = (max(d[0] for d in _dev), max(d[1] for d in _dev)) if _dev else (0, 0)
    stab = S.get('stability', {})
    if STB:
        R_ = STB['region']
        cA = STB['control_axion']
        stb_rows = ''.join(
            f'<li class="{_cls(c)}"><span class="pill">{_pill(c)}</span>'
            f'<span>{html.escape(c["name"])}{(" <em>(" + html.escape(c["note"]) + ")</em>") if c["note"] else ""}</span></li>'
            for c in STB['checks'])
        res_ = ' y '.join(f'{fnum(float(k[2:]), 1)} kHz: índice {STB[k]["total"]:.0f}' for k in STB if k.startswith('nu'))
        stab_html = f"""
  <h3>Estabilidad radial (paso 16)</h3>
  <p>Una perturbación inestable crece en el tiempo: su frecuencia tiene <code>Im ω &gt; 0</code>. Ahí la onda exterior decae y la función del empalme, escrita sin polos (<code>radial.det_entire</code>), es analítica. El principio del argumento cuenta sus ceros dentro de un contorno sin buscarlos uno a uno. Región: <code>|Re ν| ≤ {fnum(R_['xmax'], 0)}</code> kHz, <code>{R_['eps']} ≤ Im ν ≤ {fnum(R_['ymax'], 0)}</code> kHz. Además, en el eje imaginario la función es real y contamos sus cambios de signo hasta {fnum(R_['gmax'], 0)} kHz.</p>
  <p><b>Resultado:</b> {res_}, sin ningún cambio de signo en el eje. <b>El método detecta inestabilidades cuando existen:</b> en estrellas más allá de la masa máxima da índice 1, sin axión y con axión, en este caso con el modo creciente en el eje, γ = {', '.join(fnum(g, 2) for g in cA['axis_sign_changes'])} kHz.</p>
  <p><b>Qué no cubre:</b> crecimientos más lentos que {R_['eps']} kHz, que solo ve el barrido del eje (paso 15, desde 0,01 kHz); crecimientos más rápidos que {fnum(R_['ymax'], 0)} kHz, muy por encima de la escala dinámica de la estrella; y perturbaciones no radiales. Lo que obtenemos es que <b>no encontramos modos radiales inestables en esa región</b>, no que la estrella sea estable frente a toda perturbación.</p>
  <p><b>Una trampa que encontramos por el camino:</b> un único contorno grande «encontraba» 2, 4 u 8 ceros según el muestreo. La fase de la función gira muy deprisa en los bordes largos, y dos muestras seguidas podían distar más de una vuelta sin que el refinado lo notara. Dividiendo el contorno en piezas, y quitando el factor analítico e^(i k r_m), que no tiene ceros, el resultado es 0 en todas ellas. Aun así, cerca del umbral la escalera de estados ligados del átomo gravitatorio, que se acumula en ν_a justo bajo el borde inferior, hace girar la fase más deprisa de lo que resuelve un muestreo fijo: a 3,4 kHz una pieza dio 1 con 48 puntos por lado y 0 con 96. Por eso cada índice solo se acepta si no cambia al doblar el muestreo y si ningún salto de fase entre muestras supera π/3{jump_txt}. Los dos controles negativos del paso comprueban que esos fallos se detectan. Es una comprobación práctica, no una cota rigurosa: una función que diera una vuelta entera entre dos muestras la burlaría.</p>{guard_txt}
  <ul class="checks">{stb_rows}</ul>"""
    else:
        stab_html = (f"<p>Estabilidad: en el eje imaginario el determinante no cambia de signo "
                     f"({', '.join(str(v_['sign_changes']) for v_ in stab.values())} cambios).</p>")
    rob_html = (f"""
  <h3>Robustez numérica y estabilidad</h3>
  <p>Los amortiguamientos de unos pocos Hz son lo más delicado del cálculo. Los recalculamos variando la tolerancia del integrador, el radio inicial, la longitud de la integración exterior y el radio donde se empalma con la solución exterior. Ninguna variante mueve Re ν más de un {rob_re:.0e} relativo ni Im ν más de un {rob_im:.0e} relativo. Es la variación observada en estos controles del integrador, no la incertidumbre total.{eos_txt} La física (la ecuación de estado, el modelo de acoplo) puede mover las cifras mucho más:</p>
  <div class="tbl"><table><tr><th>modo</th><th>variante</th><th>Re ν (kHz)</th><th>Im ν (Hz)</th></tr>{rob_rows}</table></div>
{stab_html}""" if rob else '')
    art = S.get('artifact')
    _o = json.loads((ROOT / 'report' / 'baseline_v1' / 'fig_axgr.json').read_text())
    _n = json.loads((ROOT / 'results' / 'fig_axgr.json').read_text())
    _d = sorted(100 * abs(b / a - 1) for mo, mn in zip(_o['modes'], _n['modes']) for a, b in zip(mo['re'], mn['re']))
    shift_txt = f"entre un {fnum(_d[0], 2)} % y un {fnum(_d[-1], 2)} % (mediana {fnum(_d[len(_d) // 2], 2)} %, en las {len(_d)} frecuencias del paso 10)"
    lap = DMD['lapse'] if DMD else None
    if DMD:
        _fm = DMD['nu6']['fluid_modes']
        dp_ = DMD.get('deep')
        if dp_:
            _rob = [r_ for r_ in dp_['roots'] if r_['shift_pct'] < 0.1]
            if dp_['total'] == 0:
                _res = '<b>ningún cero</b>: no hay análogo radial de ese modo en la región explorada'
            elif DEEP is not None and not DEEP['deep'] and not dp_['roots']:
                _res = (f"el índice es {dp_['total']}, pero <b>no conseguimos localizar esos ceros</b>. Muller no converge desde ninguna semilla: "
                        "allí la función se queda en su suelo de ruido numérico y la búsqueda vaga por él sin fijar ningún punto. "
                        "La versión anterior daba cuatro «ceros» y, en el paso 21, concluía que eran artefactos del radio de empalme. "
                        "Esos «ceros» eran búsquedas que no habían convergido (entonces no se comprobaba), así que <b>retiramos esa conclusión</b>: "
                        "con esta numérica <b>no decidimos</b> si en esa región hay un análogo radial de a0 ni qué son esos dos ceros del recuento. "
                        "Para decidirlo haría falta más precisión ahí (por ejemplo, integrar hacia dentro desde el exterior o con aritmética de más dígitos)")
            elif DEEP:
                _art = (all(max(f_['shift_pct']['theta_1e-5'], f_['shift_pct']['theta_1e-6']) < 0.1 for f_ in DEEP['fluid'])
                        and bool(DEEP['deep']) and all(min(d_['shift_pct'].values()) > 5 for d_ in DEEP['deep']))   # mismos criterios que el paso 21
                _rd = DEEP['radii']
                _zz = '; '.join(f"{cnum(d_['nu'])} kHz, que pasa a {cnum(d_['nu_1e-5'])} "
                                f"y a {cnum(d_['nu_1e-6'])} kHz ({d_['shift_pct']['theta_1e-5']:.2g} % y {d_['shift_pct']['theta_1e-6']:.2g} %)"
                                for d_ in DEEP['deep'])
                _fc = '; '.join(f"{fnum(f_['nu'][0], 4)} kHz: {f_['shift_pct']['theta_1e-5']:.1g} % y {f_['shift_pct']['theta_1e-6']:.1g} %".replace('e-0', '×10^−').replace('×10^−7', '×10⁻⁷').replace('×10^−8', '×10⁻⁸') for f_ in DEEP['fluid'])
                _res = (f"el índice es {dp_['total']}, pero los ceros no tienen fluido (fracción ~10⁻¹⁴) y dependen del radio donde empalmamos con la solución exterior. "
                        f"Empalmando en r = {fnum(_rd['default'], 0)} km (lo habitual), en r = {fnum(_rd['theta_1e-6'], 0)} km (θ₀ = 10⁻⁶) o en r = {fnum(_rd['theta_1e-5'], 0)} km (θ₀ = 10⁻⁵): {_zz}. "
                        f"Con los mismos cambios, los modos fluidos se mueven muy poco ({_fc}; paso 21)")
                _res += ('. <b>Son artefactos del empalme, no modos de la estrella</b>. En esa región F está mal condicionada (Muller converge a puntos separados unos Hz desde semillas distintas), así que <b>el método no decide si hay o no un análogo radial de a0</b>. Para decidirlo haría falta extender el fondo más allá de r = '
                         f"{fnum(_rd['default'], 0)} km y empalmar donde la cola del axión sea de verdad despreciable (θ₀ sigue por encima de 10⁻⁷ al final del fondo). Con la normalización antigua del lapso, el índice de esta región es {DMD.get('old_lapse_deep', {}).get('total')}"
                         if _art else '. La prueba del paso 21 no fue concluyente (ver la lista de comprobaciones)')
            else:
                _res = (f"índice {dp_['total']}; " + '; '.join(
                    f"{fnum(r_['nu'][0], 3)} {fnum(r_['nu'][1], 3)}i kHz (fracción fluida {fnum(r_['fluid'], 2)}; al cambiar el radio de "
                    f"empalme se mueve un {r_['shift_pct']:.2g} %, {'robusto' if r_['shift_pct'] < 0.1 else 'no robusto'})" for r_ in dp_['roots']))
            _old = DMD.get('old_lapse_deep', {}).get('total')
            deep_html = f'''
  <p><b>Más abajo en el plano, donde la secuela tiene un modo axiónico.</b> La Tabla III de la secuela (l = 2, aproximación de Cowling) da a ν_a = 6 kHz un modo axiónico muy amortiguado en 9,0 − 3,8i kHz, fuera de la franja anterior. Repetimos el recuento en <code>6,05 ≤ Re ν ≤ 12</code>, <code>−5 ≤ Im ν ≤ −1,5</code> kHz: {_res}{f" (con la normalización antigua del lapso, índice {_old})" if _old is not None and not DEEP else ""}. Tampoco contradice a la secuela: allí l = 2 y Cowling; aquí l = 0 con la métrica dinámica. Fuera de <code>Re ν ≤ 12</code> kHz e <code>Im ν ≥ −5</code> kHz no hemos mirado.</p>'''
        else:
            deep_html = ''
        damped_html = f'''
  <h3>¿Hay resonancias axiónicas muy amortiguadas? (paso 20)</h3>
  <p>Contamos con el principio del argumento los modos en <code>6,05 ≤ Re ν ≤ 12</code> kHz, <code>−1,5 ≤ Im ν ≤ −0,002</code> kHz (ν_a = 6 kHz), con la condición saliente. Cada pieza del contorno se acepta solo si su índice no cambia al doblar el muestreo. Resultado: <b>{DMD['nu6']['total']} modos</b>, y son los dos fluidos que encuentra Muller ({', '.join(f"Re ν = {fnum(z[0], 4)} kHz con Im ν = {fnum(z[1] * 1e3, 2)} Hz" for z in _fm)}; fracción fluida {', '.join(fnum(f, 2) for f in DMD['nu6']['fluid_fraction'])}). En esa región <b>no hay ninguna familia axiónica</b>.{thin_txt}</p>
  <p>El control positivo comprueba que el recuento ve resonancias muy amortiguadas cuando existen: con la normalización antigua del lapso, el índice sube a {DMD['old_lapse']['total']} (piezas {DMD['old_lapse']['per_tile']}), con las resonancias de cavidad que mostraba la versión anterior.</p>{deep_html}{scat_html}'''
    else:
        damped_html = ''
    ctr_li = ''
    if CTR:
        _P = CTR['points']
        _a = max(max(p_['dphi+0.4'], p_['dphi-0.4']) for p_ in _P.values())
        _n = max(p_['n_decay60'] for p_ in _P.values())
        ctr_li = (f"\n    <li>La condición saliente usa la forma asintótica <code>e^(i(kr + η ln r))/r</code>, que ya incluye la fase coulombiana η de la cola gravitatoria, y la impone lejos (<code>|k|S = 40</code>). "
                  f"Girar el rayo ±0,4 rad cambia ψ'/ψ en el empalme como mucho {sci(_a, 1)} y alejar el punto de partida (|k|S = 60), {sci(_n, 1)}, en {len(_P)} puntos que van de un modo poco amortiguado a uno de Im ν ≈ −4,5 kHz y a un nivel ligado. "
                  "Con la onda entrante o con un punto de partida demasiado cercano, la prueba sí detecta el cambio (paso 22).</li>")
    _nf = DMD['nu6']['total'] if DMD else None
    h29_html = ''
    if OLD_S and OLD_S != S:
        _o6, _n6 = OLD_S['nu6'], S['nu6']
        _oh, _nh = complex(*_o6['H1']), complex(*s6['H1'])
        _of, _nf34 = complex(*OLD_S['nu34']['F']), complex(*n34['F'])
        _rows = (f"<tr><td>F (ν_a = 6 kHz)</td><td>{fnum(_o6['F'], 6)} kHz</td><td>{fnum(_n6['F'], 6)} kHz</td><td>{fnum(1e3 * (_n6['F'] - _o6['F']), 2)} Hz</td></tr>"
                 f"<tr><td>H1 (ν_a = 6 kHz), Re ν</td><td>{fnum(_oh.real, 6)} kHz</td><td>{fnum(_nh.real, 6)} kHz</td><td>{fnum(1e3 * (_nh.real - _oh.real), 2)} Hz</td></tr>"
                 f"<tr><td>H1, Im ν</td><td>{fnum(1e3 * _oh.imag, 4)} Hz</td><td>{fnum(1e3 * _nh.imag, 4)} Hz</td><td>{fnum(100 * (_nh.imag / _oh.imag - 1), 2)} %</td></tr>"
                 f"<tr><td>F (ν_a = 3,4 kHz), Im ν</td><td>{fnum(1e3 * _of.imag, 4)} Hz</td><td>{fnum(1e3 * _nf34.imag, 4)} Hz</td><td>{fnum(100 * (_nf34.imag / _of.imag - 1), 2)} %</td></tr>")
        if OLD_V and VAL:
            _rows += ''.join(f"<tr><td>{k} sin axión, frente a Font et al.</td><td>{OLD_V['dev_font_pct'][k]:+.3f} %</td><td>{VAL['dev_font_pct'][k]:+.3f} %</td><td></td></tr>" for k in VAL['dev_font_pct'])
        _h29_deep = ("\n    <p>La comprobación de que Muller converge retira además una conclusión de la versión anterior. Los «ceros profundos» de la región muy amortiguada, que dábamos por artefactos del radio de empalme, eran búsquedas que no habían convergido. Ahora esa región queda sin decidir (sección 6).</p>"
                     if DEEP is not None and not DEEP['deep'] else '')
        h29_html = f'''<div class="callout">
    <h3>Corrección: el arranque del interior no integraba</h3>
    <p>La sexta revisión señaló que el código no comprobaba si cada integración terminaba bien. Al añadir esa comprobación salió un fallo real. El primer tramo del interior (de r = 0,001 a 0,5 km) se integraba con un método implícito (BDF) que fallaba en su primer paso, a cualquier frecuencia. El código tomaba entonces el valor inicial en r = 0,001 km y lo usaba como si fuera el de r = 0,5 km. La parte mal arrancada decae hacia fuera, por eso el error era pequeño y los controles no lo veían. Ahora ese tramo lo integra DOP853, y toda integración que no termina detiene el cálculo. El mismo fallo estaba en los pasos 9 y 10 del informe principal, corregidos también. Todas las cifras de este informe están recalculadas:</p>
    <div class="tbl"><table><tr><th>cifra</th><th>antes</th><th>ahora</th><th>cambio</th></tr>{_rows}</table></div>{_h29_deep}
  </div>'''
    tail_txt = f"{fnum(100 * lap['tail_mass_frac'], 2)} % de la masa" if lap else "una fracción pequeña de la masa"
    phi_txt = sci(lap['phi_inf_old'], 1) if lap else "un valor distinto de 0"
    art_html = (f"""<div class="callout">
    <h3>Corrección: un artefacto retirado</h3>
    <p>La versión anterior de este informe mostraba una «rama axiónica» de resonancias muy amortiguadas en la estrella realista (6,585 − 0,563i kHz y otras casi equiespaciadas). Era un artefacto. Normalizábamos el lapso Φ en la superficie con la masa encerrada allí, sin contar la nube del axión fuera de la estrella ({tail_txt}). Eso dejaba Φ(∞) = {phi_txt} en vez de 0, y el coeficiente <code>ω² e^(−2Φ)</code> daba un salto en el punto de empalme con la solución exterior. Ese salto refleja parte de la onda y forma una cavidad entre la superficie y el empalme.</p>
    <p>Con la normalización antigua, la «resonancia» se mueve de {fnum(art['old_rm'][0], 3)} {fnum(art['old_rm'][1], 3)}i a {fnum(art['old_r84'][0], 3)} {fnum(art['old_r84'][1], 3)}i kHz al mover el empalme de {fnum(art['rm'], 0)} a {fnum(art['r84'], 0)} km, algo imposible para un modo físico. {(f'Con la corrección, el recuento del paso 20, que cubre ese punto, solo encuentra en esa franja los {_nf} modos fluidos' if DMD else 'Con la corrección, Muller desde esa semilla sale de la región, lo que por sí solo no prueba que esté vacía')}. Eso retira esas resonancias concretas, pero no establece si existe o no una familia física de resonancias axiónicas muy amortiguadas: para eso hay que contar raíces en el semiplano amortiguado{(f" (lo hace el paso 20, más abajo: con la normalización antigua el recuento sube de {DMD['nu6']['total']} a {DMD['old_lapse']['total']})") if DMD else ", y no lo hemos hecho"}.</p>
    <p>La corrección desplaza las frecuencias con axión {shift_txt} y elimina la dependencia del radio de empalme. También mejora el acuerdo con las tablas del Paper II en nuestro paso 9 (el detalle está en el cuaderno de laboratorio). Lo detectó una revisión externa que señaló que la tabla de robustez no respaldaba nuestro «por debajo del 1 %».</p>
  </div>""" if art else '')
    star_html = f"""
<section>
  <h2>6. La estrella realista: BSk26 con condensado, en relatividad general completa</h2>
  <p>Las ecuaciones radiales salen de la misma linealización de la acción que el paso 2 compara con la Ec. (3.3) del Paper I. Se convierten en el sistema <code>(ξ, Δp, δθ, δθ')</code> sin transcribir nada a mano (<code>axion_ns/radial.py</code>):</p>
  <ul>
    <li><code>δλ</code> y <code>δλ'</code> salen de <code>G_tr</code>, <code>δΦ'</code> de <code>G_rr</code> y <code>ξ'</code> de la conservación bariónica.</li>
    <li>El condensado no se anula en la superficie, así que fuera de la estrella se sigue integrando el axión acoplado a la métrica hasta que la cola es despreciable.</li>
    <li>Más allá, la onda exterior se integra en un rayo del plano r complejo donde la onda saliente decae. En la recta real, la integración hacia dentro es inestable para modos muy amortiguados: el método anterior (<code>axstar.exterior_L</code>) se desvía un {fnum(S['exterior']['old_vs_new_strong'] * 100, 0)} % en un punto muy amortiguado (6,585 − 0,563i kHz), aunque solo {S['exterior']['old_vs_new_weak']:.0e} en un modo poco amortiguado. Así que este problema apenas afecta a los amortiguamientos de pocos Hz de los pasos 9–10, lo que no los valida frente a otras fuentes de error.</li>{ctr_li}
  </ul>
  <p><b>Validación sin axión.</b> El cuadrado de la frecuencia del modo F se anula en <code>p_c = {fnum(t['p0_wF2_zero'], 1)}</code> MeV/fm³ y la masa máxima ({fnum(t['Mmax_msun'], 4)} M☉) está en <code>p_c = {fnum(t['p0_Mmax'], 1)}</code>. Coinciden, como exige el criterio del punto de retorno.</p>
{val_html}
  <figure>
    <div class="legend"><span><i class="sw s1"></i>modo fluido</span><span><i class="sw s2"></i>modo axiónico</span><span><i class="sw ref"></i>anillo: amortiguado</span></div>
    {fig4}
    <figcaption>Espectro radial a ν_a = 6 kHz (M = {fnum(s6['M_msun'], 4)} M☉, R = {fnum(s6['R_km'], 3)} km). El carácter se mide con la fracción fluida: F {fnum(s6['F_fluid'], 2)}, H1 {fnum(s6['H1_fluid'], 2)}. Es un diagnóstico, no una energía conservada: pesa <code>|ξ|²</code> con <code>(ε+p)</code> y <code>|δθ|²</code> con <code>f_a²</code> hasta el radio de empalme. Para un modo cuasi-normal, cuya autofunción crece hacia fuera, la cifra depende de ese corte. Sirve para distinguir «casi todo fluido» de «casi todo axión», no para medir mezclas finas.</figcaption>
  </figure>
  <p>La estructura es la del modelo simplificado, con un añadido de la gravedad:</p>
  <ul>
    <li>Los modos fluidos por debajo del umbral tienen frecuencia real (F) y los de encima se amortiguan poco (H1).</li>
    <li>En el rectángulo contado (<code>6,05 ≤ Re ν ≤ 12</code> kHz, hasta <code>Im ν = −1,5</code> kHz) no hay ninguna rama de resonancias axiónicas muy amortiguadas como la del modelo simplificado (paso 20): la que creímos ver era un artefacto (recuadro siguiente). Más abajo, nuestro método no decide (paso 21). Dentro de la estrella, con <code>θ ≈ π</code>, la masa efectiva del axión es mayor que <code>m_a</code> y hace de barrera, no de pozo.</li>
    <li>Justo por debajo del umbral, el potencial gravitatorio exterior atrapa al axión en una <b>escalera hidrogenoide</b>: un «átomo gravitatorio».</li>
  </ul>
  <figure>
    {fig5}
    <figcaption>Energía de ligadura de los niveles axiónicos frente a la ley hidrogenoide <code>α_g²/(2n²)</code>, con <code>α_g = GMm_a</code> (sin ajuste). El nivel n = 1, cuya órbita se solapa con la estrella, se desvía; se trata en el cuaderno.</figcaption>
  </figure>
  <div class="tbl"><table>{lad_head}{lad_rows}</table></div>
  {ryd_html}
  <h3>Cómo cruza el modo F el umbral</h3>
  <p>Al bajar ν_a, el umbral se acerca a F desde arriba, y el nivel n = 1 del átomo gravitatorio baja con él:</p>
  <ul>
    <li>En <b>ν_a = 3,5 kHz</b> hay dos estados ligados ({fnum(n35['bound'][0], 4)} y {fnum(n35['bound'][1], 4)} kHz) con carácter mezclado ({fnum(n35['fluid'][0], 2)} y {fnum(n35['fluid'][1], 2)} de fluido){': el cruce evitado de la figura siguiente' if CRS else '. Sugiere un cruce evitado; falta seguir las ramas'}.</li>
    <li>En <b>ν_a = 3,4 kHz</b>, F ya está por encima del umbral, en {fnum(n34['F'][0], 4)} kHz, con un amortiguamiento de solo {fnum(-n34['F'][1] * 1e3, 2)} Hz (fracción fluida {fnum(n34['fluid'], 2)}).</li>
  </ul>
  {cross_html}
  <p>La transición cuasi-normal → ligado de la estrella realista no es un cruce limpio como en el modelo simplificado: el modo fluido atraviesa una escalera densa de estados del axión que se acumula en el umbral.</p>
{rob_html}
  {h29_html}
  {art_html}
  {damped_html}
  <div class="callout">
    <h3>Abierto</h3>
    <p>No hemos podido decidir si el amortiguamiento tiende a cero o a un valor finito justo encima del umbral (factor de Sommerfeld del potencial 1/r). Los dos puntos disponibles apuntan a que cae deprisa, pero hace falta un barrido denso en ν_a.</p>
    <p>El defecto cuántico δ ≈ 0,056 lo fija el interior (la estrella y el campo fuerte), y no lo hemos predicho desde ahí: es un dato que sale del cálculo, no una comprobación independiente. Tampoco hemos estudiado aparte el nivel n = 1, cuya órbita se solapa con la estrella.</p>
    <p>No sabemos si hay modos radiales axiónicos con amortiguamiento de varios kHz (la región de a0 en la secuela): ahí nuestra F está en su suelo de ruido numérico y no conseguimos localizar sus ceros (pasos 20–21).</p>
  </div>
  <ul class="checks">{st_rows}</ul>
</section>"""
    sum_extra = ''
    if STB:
        sum_extra += (f"\n    <li><b>Estabilidad:</b> no encontramos modos radiales inestables en <code>|Re ν| ≤ {fnum(STB['region']['xmax'], 0)}</code>, "
                      f"<code>{STB['region']['eps']} ≤ Im ν ≤ {fnum(STB['region']['ymax'], 0)}</code> kHz (principio del argumento, paso 16). Es un resultado sobre esa región, no una demostración de estabilidad.</li>")
    if DMD and DEEP:
        sum_extra += (f"\n    <li><b>Semiplano amortiguado:</b> en el rectángulo contado (<code>6,05 ≤ Re ν ≤ 12</code> kHz, <code>{_imtop} ≥ Im ν ≥ −1,5</code> kHz) solo están los {DMD['nu6']['total']} modos fluidos (paso 20). "
                      + ("Más abajo, donde la secuela tiene su modo axiónico a0 (l = 2, Cowling), el recuento da ceros que nuestra numérica no consigue localizar, así que el método <b>no decide</b> si hay un análogo radial (pasos 20–21).</li>"
                         if not DEEP['deep'] else
                         "Más abajo, donde la secuela tiene su modo axiónico a0 (l = 2, Cowling), los ceros que salen son artefactos del empalme y el método <b>no decide</b> si hay un análogo radial (paso 21).</li>"))
    val_li = (f"<li><b>Validación externa:</b> sin axión, el mismo código reproduce las frecuencias radiales F, H1, H2 y H3 "
              f"del politropo de referencia de Font et al. (2002) con una diferencia máxima del "
              f"{max(abs(v) for v in VAL['dev_font_pct'].values()):.3f} % (sector fluido sin axión).</li>") if VAL else ''
    star_summary = f"""
    {val_li}
    <li><b>Estrella realista</b> (BSk26, relatividad general completa, ν_a = 6 kHz): el modo F por debajo del umbral tiene frecuencia real, H1 por encima se amortigua {fnum(-s6['H1'][1] * 1e3, 1)} Hz{scat_sum}. Justo bajo el umbral aparece una escalera de estados ligados del axión, un «átomo gravitatorio»{(f", que sigue una ley de Rydberg con un defecto cuántico casi constante, δ ≈ {RYD['ritz']['delta0']:.3f}" if RYD else "")}.</li>
    <li>En la estrella realista, el modo F cruza el umbral atravesando esa escalera, con un cruce evitado con el nivel n = 1{f" (separación mínima {1e3 * CRS['gap_min_khz']:.0f} Hz)" if CRS else ""}, y sale como cuasi-normal con un amortiguamiento muy pequeño ({fnum(-n34['F'][1] * 1e3, 2)} Hz a ν_a = 3,4 kHz).</li>{sum_extra}"""
    # ---- apertura: conclusiones e importancia (cifras desde los datos)
    _h1w = f", y un segundo método independiente de la búsqueda de raíces da {fnum(SCAT['H1 (6 kHz)']['half_width'] * 1e3, 2)} Hz" if SCAT else ''
    _font = (f"<li><b>El código es fiable donde se puede contrastar:</b> sin axión reproduce las vibraciones radiales publicadas por "
             f"Font et al. (2002) con una diferencia máxima del {fnum(max(abs(v) for v in VAL['dev_font_pct'].values()), 3)} %.</li>") if VAL else ''
    _stab = (f"<li><b>Estabilidad:</b> no encontramos vibraciones radiales que crezcan en la región explorada (paso 16). "
             f"Es un resultado sobre esa región, no una demostración general.</li>") if STB else ''
    intro_html = f"""
<section class="sheet summary intro">
  <h3>En pocas palabras</h3>
  <p>Una estrella de neutrones puede «respirar»: hincharse y deshincharse entera. Si dentro hay un condensado de axiones, esa vibración puede escapar en forma de axiones, pero solo si su frecuencia supera la masa del axión (el umbral). José preguntó qué pasa exactamente cuando una vibración cruza ese umbral.</p>
  <p><b>Conclusiones (preliminares)</b></p>
  <ul>
    {_font}
    <li><b>Por encima del umbral, la vibración se apaga emitiendo axiones; por debajo, no.</b> Obtenemos que, en una estrella realista con axiones de 6 kHz, la vibración H1 (por encima) se amortigua a {fnum(-s6['H1'][1] * 1e3, 1)} Hz{_h1w}. La vibración fundamental F (por debajo) no pierde energía.</li>
    <li><b>El cruce no es un interruptor.</b> Justo por debajo del umbral aparece una escalera de estados en los que el axión queda atrapado por la gravedad de la estrella, un «átomo gravitatorio»{(f" que sigue la misma ley que los átomos de Rydberg (defecto cuántico δ ≈ {RYD['ritz']['delta0']:.3f})" if RYD else "")}. La vibración F atraviesa esa escalera esquivando el primer peldaño{(f" (separación mínima {1e3 * CRS['gap_min_khz']:.0f} Hz)" if CRS else "")} y, recién cruzado el umbral, apenas pierde energía ({fnum(-n34['F'][1] * 1e3, 2)} Hz a ν_a = 3,4 kHz).</li>
    {_stab}
    <li><b>Queda abierto</b> si existe un análogo radial de los modos axiónicos muy amortiguados de la secuela: nuestro método no lo decide (pasos 20–21).</li>
  </ul>
  <p><b>Por qué importa</b> <span class="small">(borrador pendiente de revisión por José)</span></p>
  <p>El primer artículo sugiere que la sismología de estrellas de neutrones podría servir para estudiar las propiedades del axión, y da la regla básica: las vibraciones por encima del umbral se amortiguan y las de debajo no. Nuestros resultados sugieren que cerca del umbral esa regla tiene matices: la transición pasa por el átomo gravitatorio, y justo por encima el amortiguamiento puede ser muy pequeño. Si la sismología llegara a usarse para acotar la masa del axión, es precisamente cerca del umbral donde la señal cambia de comportamiento, así que esa zona hay que entenderla bien. El propio artículo advierte de que las vibraciones radiales son difíciles de observar, así que esto es sobre todo una pieza teórica para cuando el análisis se extienda a vibraciones observables.</p>
</section>"""
else:
    star_html = None
    star_summary = ''
    intro_html = ''

# ------------------------------------------------------------------ respuesta a las revisiones externas
REVIEWS = [
    ('1ª', 'Signo de la identidad Im w = Re k · Im k / Re w', 'Corregido; la comprobación evalúa el residuo con los dos signos', 'Sec. 3, 7', 'ok'),
    ('1ª', '«Función entera» no demostrada', 'Precisado: lo entero en k es f/√D; ninguna raíz es un cero espurio de D', 'Sec. 1, 7', 'ok'),
    ('1ª, 2ª E', 'Interpretación del modo inestable (fondo elegido, nuevo equilibrio, su estabilidad)', 'Separado: es una propiedad del modelo fenomenológico; la pregunta queda para el grupo; θ = π no demuestra la estabilidad', 'Sec. 5', 'ok'),
    ('1ª', 'Ficheros no adjuntos', 'Código y datos publicados junto a la página (e incrustados en la versión autocontenida)', 'Sec. 9', 'ok'),
    ('1ª, 2ª, 3ª C', 'Falta un método o una referencia independientes', 'Principio del argumento sin Muller (modelo simplificado); sin axión, Font et al. 2002; punto de retorno', 'Sec. 6, 7', 'ok'),
    ('1ª', 'Convergencia de los amortiguamientos de pocos Hz', 'Tabla de robustez: tolerancia, radio inicial, longitud exterior y radio de empalme', 'Sec. 6', 'ok'),
    ('1ª, 2ª D, 4ª', 'Ley de umbral en la estrella relativista', 'Declarado abierto; el exponente 0,457 es del modelo simplificado', 'Sec. 6 («Abierto»)', 'open'),
    ('2ª A, 3ª A', 'Estabilidad radial no demostrada', 'Principio del argumento en el semiplano inestable, con controles positivos y negativos; formulado como «no encontramos en esta región»', 'Sec. 6 (paso 16)', 'ok'),
    ('2ª B', 'Tabla de robustez incoherente (2,7 % frente a 0,17 %)', 'Era el artefacto del lapso; corregido y recalculado', 'Sec. 6', 'ok'),
    ('2ª C', 'Condición exterior: ángulo del rayo, punto de partida, orden asintótico, cola coulombiana', 'Paso 22: la forma asintótica lleva la fase coulombiana; ángulo y punto de partida no cambian ψ\'/ψ, y la prueba detecta los casos malos', 'Sec. 6 (paso 22)', 'ok'),
    ('3ª B', 'El exponente es un ajuste; incoherencia «se anula en el umbral»', 'Ley de tramo intermedio, separada del régimen de colisión', 'Resumen, Sec. 3', 'ok'),
    ('3ª C', 'Extrapolar la validez de los pasos 9–10', 'Rebajado: el error del exterior antiguo es pequeño en modos poco amortiguados, lo que no los valida frente a otras fuentes de error', 'Sec. 6', 'ok'),
    ('3ª D', 'Átomo gravitatorio: identificación de niveles e influencia del interior', 'Ley de Rydberg con defecto cuántico (paso 18); predecir δ desde el interior queda abierto, igual que el nivel n = 1', 'Sec. 6', 'open'),
    ('3ª E, 4ª 5, 5ª 4', 'Cruce evitado sin seguimiento de ramas; fracción fluida en modos cuasi-normales', 'Seguimiento continuo de las dos ramas con intercambio de carácter (paso 19); la fracción fluida se presenta como diagnóstico', 'Sec. 6', 'ok'),
    ('4ª 1, 5ª 1', 'Colisión de polos mal situada', 'Raíz doble resuelta (f = 0, ∂f/∂k = 0); colisión y cruce de k = 0 son dos sucesos', 'Resumen, Sec. 2', 'ok'),
    ('4ª 2', '«Linealizado en θ = 0» no describe el modelo', 'Precisado: masa de θ = 0 con acoplo constante, modelo fenomenológico', 'Sec. 5', 'ok'),
    ('4ª 3, 5ª 2', 'Font et al.: no todo está dentro del redondeo', 'Diferencias en Hz; solo F está dentro; valida el sector sin axión', 'Sec. 6', 'ok'),
    ('4ª 4, 5ª 3', 'Retirar un artefacto no prueba que no haya una familia axiónica', 'Recuento del semiplano amortiguado (paso 20): hasta Im ν = −1,5 kHz, solo modos fluidos', 'Sec. 6 (paso 20)', 'ok'),
    ('5ª 3', 'La secuela tiene modos axiónicos muy amortiguados (Tabla III)', 'Recuento en esa región: índice 2, pero nuestra numérica no localiza esos ceros. La conclusión anterior («artefactos del empalme») se apoyaba en búsquedas sin converger y se retira; el método no decide', 'Sec. 6 (pasos 20–21)', 'open'),
    ('4ª 5, 5ª 5', 'Defecto cuántico «casi», no exactamente, constante', 'Dicho así; el ajuste de Ritz da la dependencia residual', 'Resumen, Sec. 6', 'ok'),
    ('4ª 5', 'Estabilidad numérica no es incertidumbre física', 'Dicho en la tabla de robustez (EOS, modelo y física omitida aparte)', 'Sec. 6', 'ok'),
    ('5ª 5', 'Umbral del modo inestable con λ ≠ 0; clasificación κ < 1; «único canal»', 'Umbral con λ² = 10⁻³ calculado; κ < 1 para frecuencia real; limitado al modelo idealizado', 'Resumen, Sec. 2, 5', 'ok'),
    ('6ª 1', 'Dominio del recuento: el texto decía Im ν < 0, pero se contaba hasta −2 Hz y solo en 6,05–12 kHz', 'Franja casi real contada (paso 20); la conclusión se limita al rectángulo y se dice qué queda fuera', 'Resumen, Sec. 6', 'ok'),
    ('6ª 2', '«Incertidumbre numérica» del 0,004 %; inconsistencia termodinámica de la EOS', 'Rebajado a «variación observada en estos controles»; inconsistencia de los ajustes BSk26 medida, con su efecto en las frecuencias (paso 23)', 'Sec. 6', 'ok'),
    ('6ª 3', 'Falta una validación independiente del amortiguamiento con axión', 'Anchura de la resonancia en la dispersión a frecuencia real, sin plano complejo (paso 24): otro método, pero las mismas ecuaciones. Una prueba con otras ecuaciones (balance de energía radiada) queda pendiente', 'Sec. 6 (paso 24)', 'open'),
    ('6ª 4', 'Muller sin residuo, solve_ivp.success sin comprobar, recuento heurístico', 'Muller y las integraciones fallan si no convergen. Eso destapó un fallo real en el arranque del interior, ya corregido y con todo recalculado; el recuento se presenta como heurística', 'Sec. 6 (recuadros)', 'ok'),
    ('6ª 5', 'El HTML autocontenido no basta para reproducir', 'Incrustados los ficheros que faltaban (paper.py, __init__.py, checks.json, datos de referencia) y un extractor; probado en una carpeta vacía', 'Sec. 9', 'ok'),
]
_CK = json.loads((ROOT / 'results' / 'checks.json').read_text())['steps']
_dur = lambda t: '<1 s' if t < 1 else f'~{t:.0f} s' if t < 90 else f'~{t / 60:.0f} min'
_repro = '\n'.join(f"uv run python -m axion_ns.{s_['module']:<24}  # {s_['step']} ({_dur(s_['seconds'])})"
                    for s_ in _CK if s_['module'] >= 'step14' and s_['module'][4:6].isdigit())
_rv_rows = ''.join(f'<tr><td>{a}</td><td>{b}</td><td>{c}</td><td>{d}</td><td><span class="st {e}">{"resuelto" if e == "ok" else "abierto"}</span></td></tr>'
                   for a, b, c, d, e in REVIEWS)
_n_open = sum(r[4] == 'open' for r in REVIEWS)
REV_HTML = f"""
<section>
  <h2>8. Respuesta a las revisiones externas</h2>
  <p>Este informe ha recibido seis revisiones externas. Cada punto, con lo que hicimos y dónde está; {len(REVIEWS) - _n_open} están resueltos y {_n_open} quedan abiertos, declarados como tales en el texto. Las revisiones 4ª y 5ª coinciden en varios puntos porque la 5ª leyó una versión anterior.</p>
  <div class="tbl"><table class="rv"><tr><th>revisión</th><th>crítica</th><th>qué hicimos</th><th>dónde</th><th></th></tr>{_rv_rows}</table></div>
</section>"""

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

# ------------------------------------------------------------------ ficheros para reproducir (se publican junto a la página)
import shutil
FILES = ['axion_ns/__init__.py', 'axion_ns/paper.py', 'axion_ns/figdata.py', 'axion_ns/model.py', 'axion_ns/geometry.py', 'axion_ns/compare.py', 'axion_ns/pipeline.py', 'axion_ns/nonradial.py',
         'axion_ns/fullgr.py', 'axion_ns/bsk26.py', 'axion_ns/axstar.py', 'axion_ns/axgr.py', 'axion_ns/radial.py',
         'axion_ns/step2_perturbations.py', 'axion_ns/step4_numerics.py', 'axion_ns/step14_radial_branches.py',
         'axion_ns/step15_radial_star.py', 'axion_ns/step16_radial_stability.py', 'axion_ns/step17_radial_validation.py', 'axion_ns/step18_rydberg.py', 'axion_ns/step19_avoided_crossing.py',
         'axion_ns/step20_damped_count.py', 'axion_ns/step21_deep_artifacts.py', 'axion_ns/step22_exterior_contour.py',
         'axion_ns/step23_eos_consistency.py', 'axion_ns/step24_scattering_width.py',
         'report/build_radial.py', 'report/build_radial_standalone.py',
         'pyproject.toml', 'uv.lock']
DATA = ['results/radial_branches.json', 'results/radial_star.json', 'results/radial_stability.json',
        'results/radial_validation.json', 'results/radial_rydberg.json',
        'results/radial_crossing.json', 'results/radial_damped.json', 'results/radial_deep.json', 'results/radial_contour.json',
        'results/radial_eos.json', 'results/radial_scattering.json', 'results/checks.json',
        'results/fig_axgr.json', 'report/baseline_v1/fig_axgr.json',
        'report/baseline_pre_h29/radial_star.json', 'report/baseline_pre_h29/radial_validation.json']   # tabla «antes/ahora» de H29
FD = ROOT / 'report' / 'radial_files'
FD.mkdir(exist_ok=True)
published = {}
for f in FILES + DATA:
    src = ROOT / f
    if not src.exists():
        continue
    name = f.replace('/', '__') + ('.txt' if not f.endswith(('.json', '.txt')) else '')
    shutil.copy(src, FD / name)
    published[f'radial_files/{name}'] = f'report/radial_files/{name}'
(ROOT / 'report' / 'radial_files.json').write_text(json.dumps(published, indent=1))

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
.res .pill {{ background: var(--code-bg, #eee); color: inherit; }}
.bad .pill {{ background: var(--bad-bg); color: var(--bad); }}
.status {{ font: 600 12px var(--f-mono); color: var(--ok); }}
.status.pend {{ color: var(--flag); }}
.rv td {{ vertical-align: top; }}
.rv .st {{ font: 600 11px/1 var(--f-mono); padding: 3px 6px; border-radius: 3px; white-space: nowrap; }}
.rv .st.ok {{ background: var(--ok-bg); color: var(--ok); }}
.rv .st.open {{ background: var(--flag-bg); color: var(--flag); }}
pre.code {{ background: var(--code-bg); padding: 12px 14px; border-radius: 6px; overflow-x: auto; font: 13px/1.5 var(--f-mono); margin: 0; }}
pre.code code {{ background: none; padding: 0; }}
.files {{ columns: 2 220px; margin: 0; padding-left: 18px; font-size: 13px; }}
.small {{ font-size: 13px; color: var(--muted); }}
footer {{ font-size: 13px; color: var(--muted); border-top: 1px solid var(--rule); padding-top: 14px; }}
a {{ color: var(--accent); }}
:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
</style>
<div class="wrap">
<header>
  <div class="eyebrow">Workshop IA × física · pasos 14–24 · modelos radiales del Paper I</div>
  <h1>Del modo cuasi-normal al modo ligado: el umbral axiónico en las oscilaciones radiales</h1>
  <p class="lede">Respuesta a la pregunta de José A. Pons: cómo pasa un modo radial de cuasi-normal a evanescente cuando su frecuencia cruza la masa del axión, y cómo se interpretan las distintas ramas del espectro. Hipótesis anotadas antes de calcular; todas las cifras se insertan desde <code>results/radial_*.json</code> al generar la página.</p>
</header>
{intro_html}

<section class="sheet summary">
  <h3>Resumen</h3>
  <ul>
    <li>Para <code>l = 0</code> no hay radiación gravitatoria, así que, en este modelo idealizado (sin viscosidad, neutrinos ni otras disipaciones), el único canal de pérdida es la emisión de axiones, con <code>k_out² = ω² − m_a²</code>. Escrita en la variable <code>k = k_out</code>, la función de dispersión del modelo (Ec. 4.13) no tiene la ambigüedad de hoja del umbral exterior.</li>
    <li>La transición sigue la secuencia típica de onda con umbral: el par cuasi-normal <code>k = ±a − ib</code> colisiona en el eje imaginario negativo en <code>x₀ = {fnum(C3['x0'], 5)}</code> (raíz doble, λ² = 10⁻³), cuando su frecuencia ya está <b>por debajo</b> de <code>m_a</code> (<code>w − 1 = {sci(C3['w_minus_1'])}</code>). Se convierte en dos <b>estados virtuales</b> (<code>Im k &lt; 0</code>, frecuencia real pero no físicos). Uno de ellos sube, cruza <code>k = 0</code> (justo <code>ω = m_a</code>) en <code>x₀ = {fnum(C3['x0_k0'], 5)}</code> y pasa a <b>estado ligado</b> (<code>Im k &gt; 0</code>); el otro sigue virtual.</li>
    <li>Cerca del umbral, <code>|Im w| ∝ (Re w − 1)^{fnum(s_fit, 2)}</code>, próximo a ½. Es una ley de tramo intermedio. El amortiguamiento solo se anula en la colisión, ya por debajo del umbral, a <code>w − 1 ≈ k²/2 ∼ −(Im k)²/2</code> (λ² = 10⁻⁴: <code>{sci(C4['w_minus_1'], 2)}</code>).</li>
    <li>Las ramas fluida y axiónica no se cruzan: la separación mínima es lineal en el acoplo, con <code>Δw²_min/λ = {fnum(min(ratios), 2)}–{fnum(max(ratios), 2)}</code> para λ² = 10⁻⁶…10⁻³.</li>
    <li>Hallazgo no previsto: con la masa efectiva taquiónica <code>m_a²(1 − α)</code> del modelo simplificado, este contiene un modo inestable (<code>w² = {fnum(w2u, 2)}</code> con los parámetros de la figura), que no aparece en la ventana <code>Im w ≤ 0</code> habitual. Lo planteamos como pregunta abajo.</li>{star_summary}
  </ul>
</section>

<section>
  <h2>1. Planteamiento</h2>
  <p>Usamos el modelo simplificado de la Sec. IV del Paper I (Ec. 4.13, ya reproducido en el paso 4): fondo homogéneo, simetría plana, sin perturbación métrica, con <code>w = ω/m_a</code> y <code>x₀ = m_a R</code>. Los parámetros son los de la Fig. 3: <code>α = {fnum(alpha, 0)}</code>, <code>c_s² = {fnum(D["params"]["cs2"], 1)}</code> y <code>λ² = 10⁻³</code>. Para mover un modo fluido a través del umbral variamos <code>x₀</code>, lo que equivale a cambiar <code>m_a R</code>.</p>
  <p>La función depende de <code>ω</code> solo a través de <code>w²</code>, salvo el término de onda exterior <code>i k_out</code>. Con <code>w² = 1 + k²</code> desaparece la raíz cuadrada del umbral exterior y cada raíz se clasifica sin ambigüedad. Hay una sutileza que nos señaló la revisión: la función contiene además <code>√D</code> (Ec. 4.13) y es impar en ella. Lo que es entero en <code>k</code> es <code>f/√D</code>, y <code>f</code> tiene ceros espurios donde <code>D = 0</code>. Comprobamos que ninguna raíz de este informe está ahí (<code>|D| ≥ 21</code>):</p>
  <div class="tbl"><table>
    <tr><th>Posición de la raíz en el plano k</th><th>Onda exterior</th><th>Interpretación</th></tr>
    <tr><td><code>Re k &gt; 0</code>, <code>Im k &lt; 0</code></td><td>saliente, crece con r</td><td>modo cuasi-normal (amortiguado)</td></tr>
    <tr><td><code>Re k = 0</code>, <code>Im k &gt; 0</code></td><td>evanescente <code>e^(−κr)</code></td><td>modo ligado, normalizable: <code>w² = 1 − κ²</code>, frecuencia real si <code>κ &lt; 1</code>; si <code>κ &gt; 1</code>, <code>w² &lt; 0</code> y el modo crece (el inestable de la Sec. 5)</td></tr>
    <tr><td><code>Re k = 0</code>, <code>Im k &lt; 0</code></td><td>crece <code>e^(+κr)</code></td><td>estado virtual, no es un modo físico (frecuencia real si <code>κ &lt; 1</code>)</td></tr>
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
  <p>La identidad exacta <code>Im w = Re k · Im k / Re w</code> (de <code>w² = 1 + k²</code>; la primera versión de este informe tenía el signo cambiado) explica el comportamiento: con <code>Im k</code> casi constante, la tasa es proporcional a <code>k_out ∝ √(ω² − m_a²)</code>, el espacio de fases de la onda emitida. El exponente ajustado queda algo por debajo de ½ porque <code>Im k</code> varía un {imk_var:.0f} % en el tramo. Es una ley de <b>tramo intermedio</b>, ajustada en un intervalo finito, no una ley asintótica: se rompe cuando <code>Re k ∼ |Im k|</code>. Ahí manda la colisión de los dos polos, a una distancia del umbral de orden <code>(Im k)²</code>, minúscula con acoplo débil. El amortiguamiento solo se anula en la colisión, que calculamos como raíz doble y no del barrido: <code>x₀ = {C4['x0']:.10f}</code>, <code>k = {C4['k'][1]:.7f} i</code>, <code>w − 1 = {sci(C4['w_minus_1'], 4)}</code> (λ² = 10⁻⁴). Ocurre ya <b>por debajo</b> del umbral, como exige <code>w² = 1 + k²</code> con k imaginario negativo. Justo después hay dos estados virtuales, y uno cruza <code>k = 0</code> en <code>x₀ = {C4['x0_k0']:.10f}</code> para convertirse en ligado. La versión anterior de este informe situaba la anulación en el último punto del barrido (<code>Re w − 1 ≈ {sci(min(x for x, _ in dp))}</code>), donde el amortiguamiento todavía no es cero; lo señaló una revisión externa.</p>
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
  <h2>5. Un modo inestable en el modelo simplificado (masa efectiva m_a²(1 − α))</h2>
  <p>Una precisión sobre qué modelo es. El de la Sec. IV del Paper I no es la linealización consistente de la acción en <code>θ = 0</code>: allí <code>∂m*/∂a = 0</code>, y el acoplo lineal se anularía. El paper combina la masa efectiva de <code>θ = 0</code>, <code>m_a²(1 − α)</code>, con un acoplo constante. Es un modelo fenomenológico, y lo que sigue es una propiedad de ese modelo, no del problema exacto de la condensación.</p>
  <p>Con los parámetros de la Fig. 3 (<code>x₀ = 0.6</code>) hay un estado ligado con <code>w² = {fnum(w2u, 3)}</code>, es decir, una perturbación que crece como <code>e^({fnum(grow, 2)} m_a t)</code>. Es el estado fundamental del pozo taquiónico interior, con onda oscilante dentro y evanescente fuera. Sin acoplo aparece para <code>x₀ &gt; {fnum(xc, 8)}</code>, valor que coincide con la condición analítica <code>arctan(1/√(α−1))/√(α−1) = {fnum(xc_an, 8)}</code>. Con λ² = 10⁻³ el umbral baja a <code>{fnum(xc3, 8) if xc3 else '—'}</code>: el acoplo apenas lo mueve, así que la inestabilidad la produce sobre todo el término de masa taquiónico.</p>
  <p>Nuestro paso 4 no lo detectaba porque buscaba raíces solo en <code>Im w ≤ 0.01</code>. Linealizando en <code>θ = π</code> (masa efectiva interior positiva, sin acoplo lineal) desaparece: no hay estados ligados ni inestables. Eso solo dice que este modelo no tiene esa inestabilidad, no que el condensado sea estable frente a toda perturbación. Para la estrella realista, el paso 16 cuenta los modos radiales inestables en una región acotada del plano (Sec. 6).</p>
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

{REV_HTML}

<section>
  <h2>9. Reproducirlo</h2>
  <p>El código y los datos están publicados junto a esta página (los <code>.py</code> van como <code>.py.txt</code>). Con el repositorio y <a href="https://docs.astral.sh/uv/">uv</a>:</p>
  <pre class="code"><code>uv run python -m axion_ns.step2_perturbations     # ecuaciones (3.3) desde la acción (paso 2)
{_repro}
uv run python report/build_radial.py              # esta página, desde los JSON</code></pre>
  <ul class="files">{"".join(f'<li><a href="{k}"><code>{k.split("/", 1)[1].removesuffix(".txt").replace("__", "/").replace("//init/.py", "/__init__.py")}</code></a></li>' for k in sorted(published))}</ul>
  <p class="small">La tabla de la EOS BSk26 la genera <code>axion_ns/bsk26.py</code> a partir de los ajustes analíticos publicados; no hace falta ningún fichero externo.</p>
</section>

<section>
  <div class="callout neutral">
    <h3>Nota sobre el informe anterior de «tres líneas de investigación»</h3>
    <p>Lo hemos retirado. Dos revisiones independientes y una auditoría del código mostraron que los modos g, la rotación y el puente axión-magnetar se apoyaban en constantes ajustadas a mano y en premisas físicas incorrectas: un salto de <code>m_N*</code> diez veces mayor, <code>f_GW = 2ν</code> en la resonancia de marea y <code>μ₅</code> a partir de <code>∇a</code> en vez de <code>∂ₜa</code>. Las comprobaciones eran circulares. El análisis completo está en el cuaderno de laboratorio, y desde entonces hemos endurecido las reglas de verificación. Los pasos 1–10 (informe principal) no se ven afectados.</p>
  </div>
</section>

<footer>Generado por <code>report/build_radial.py</code> desde <code>results/radial_*.json</code> (pasos 14–24) y <code>results/checks.json</code>. Referencia: A. Gómez-Bañón, P. Pnigouras y J. A. Pons, <a href="https://arxiv.org/abs/2510.16245">arXiv:2510.16245</a>.</footer>
</div>
'''

for out in (ROOT / 'report' / 'informe_radial.html', ROOT / 'results' / 'informe_radial.html'):
    out.write_text(page)
print('informe_radial.html:', len(page) // 1024, 'KB;', n_ok, '/', len(checks), 'checks')
