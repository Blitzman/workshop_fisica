"""Genera results/informe_axiones.html a partir de results/checks.json y las figuras."""
import html
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / 'results'

# Etiquetas legibles (LaTeX) en el mismo orden que las comprobaciones de cada paso.
# 'obs' marca comprobaciones que pasan pero señalan algo a comentar con los autores.
PRETTY = [
    [r'Identidad \(\nabla_\mu T^{\mu}{}_{r}[a] = (\Box a + V^\prime)\,\partial_r a\)',
     r'\(G_{tr}\) del fondo estático se anula idénticamente',
     r'Ec. (2.6a), equilibrio hidrostático \(p^\prime\)',
     r'Ec. (2.6b), \(\Phi^\prime\)',
     r'Ec. (2.6c), \(M^\prime\)',
     r'Ec. (2.6d), Klein-Gordon',
     r'Sin axión se recupera el TOV estándar',
     (r'\(p(n,a)=p(n,0)\) es exacta con la Ec. (2.7); la Ec. (2.8) la da como aproximación', 'obs')],
    [r'Ec. (3.3a), Klein-Gordon perturbada',
     r'Ec. (3.3b), componente \(rr\)',
     r'Ec. (3.3c), componente \(tr\)',
     r'Ec. (3.3d), Euler',
     r'Ec. (3.3e), continuidad',
     r'Control negativo: \(4\pi\to8\pi\) en (3.3b) se detecta',
     r'Control negativo: signo del término \(2\delta\lambda\) en (3.3a) se detecta',
     r'Control negativo: término \(\delta\lambda\) omitido en (3.3e) se detecta',
     r'\(G_{tt}\) linealizada es redundante (Bianchi + conservación bariónica)',
     r'Exterior: (3.3a) se reduce a KG en Schwarzschild, Ec. (3.7)'],
    [r'Ec. (4.2a) se sigue de (3.3d) + (3.3e) con las aproximaciones (4.1)',
     r'Ec. (4.2b) se sigue de (3.3a) con las aproximaciones (4.1)',
     r'Forma adimensional (4.5) y matriz \(\mathbb{M}\), Ec. (4.6)',
     r'\(\det(\mathbb{M}-\bar k^2) = 0\) es la Ec. (4.7)',
     r'Cuadrática en \(\bar\omega^2\), Ec. (4.8), y sus raíces',
     r'Raíces \(\bar\omega^2(\bar k)\) del texto',
     r'\(\bar k_\pm^2\) y \(\Delta\), Ec. (4.12)',
     r'Casando interior y exterior con \(\delta p(\bar x_0)=0\) sale \(f(\bar\omega)\), Ec. (4.13)',
     r'Límite \(\lambda\to0\) de (4.13) da la Ec. (4.11)',
     r'\(V^{\prime\prime}+n\,m^{*\prime\prime}=m_a^2(1-\alpha)\) en \(\theta=0\): es la aproximación (4.1b)',
     (r'\(\partial_a m^* \propto \sin\theta\): se anula en \(\theta=0\) y en \(\theta=\pi\)', 'obs'),
     (r'En \(\theta=\pi\): \(m_\mathrm{eff}^2 = m_a^2(\alpha-1)/\sqrt{1-\beta} > 0\)', 'obs')],
    [r'\(\beta = 4z/(1+z)^2 = 0.877\) (paper: 0.88)',
     r'Reducción de masa \(\sigma_N[1-g(\pi f_a)] = 32.4\) MeV (paper: ~32 MeV)',
     r'\(m_a = 3\) kHz \(\leftrightarrow\) 12.4 peV (paper: 12 peV)',
     r'\(m_a R\) para 3 kHz y 12 km: 0.75 (paper: 0.12)',
     r'\(m_a^2 f_a^2\): \(\epsilon\times\)4.4–4.8 MeV/fm³ (paper: \(\epsilon\times\)13)',
     r'\(\lambda^2\) con \((\varepsilon+p)=500\) MeV/fm³, \(c_s^2=0.3\): \(\epsilon\times\)(0.029–0.032) (paper: \(\epsilon\times\)0.02)',
     r'\(n_c(\epsilon=1) = 0.40\) fm\(^{-3}\) \(\approx 2.5\,n_\mathrm{sat}\) (paper: ~2 \(n_\mathrm{sat}\))',
     r'Fig. 3 reproducida: 10 modos fluidos + 5 axiónicos en las mismas posiciones',
     (r'Variante en \(\theta=\pi\): modos fluidos sin amortiguar', 'obs'),
     r'Figura generada'],
    # ---- Paper II (arXiv:2610.37308)
    [r'Ec. (8a), momento radial',
     r'Ec. (8b), conservación bariónica',
     r'Ec. (8c), momento angular: \((\varepsilon+p)\mathcal{V}\) algebraico',
     r'Ec. (8d), Klein-Gordon interior, con \(m_\mathrm{eff}^2\) de la Ec. (9)',
     r'Control negativo: quitar el \(l/r\) de \(\delta\theta^\prime + l\,\delta\theta/r\) en (8a) se detecta',
     r'Control negativo: quitar \(\Phi^\prime\) del lado izquierdo de (8b) se detecta',
     r'Control negativo: \(2l/r\) en vez de \(2(l+1)/r\) en (8d) se detecta',
     r'Exterior: (8d) con \(n_b=0\) y \(\delta\theta=(R/r)^l\delta\theta_\mathrm{ext}\) es la Ec. (13)',
     r'Ec. (15), ajuste de derivadas logarítmicas en \(R\)',
     r'Ec. (10b): el orden \(1/r\) de (8d) impone \(\delta\theta^\prime(0)=0\)',
     r'Ec. (10a): \(\delta p(0)\) del orden \(1/r\) de (8a)',
     r'El orden \(1/r\) de (8b) da la misma condición que (10a)',
     r'Ecs. (16a,b): \(\delta p^\prime(0)=0\) y \([(\varepsilon+p)W]^\prime(0)=0\)',
     r'Ec. (16c): \(\delta\theta^{\prime\prime}(0)\), con su \(1/(2l+3)\)',
     r'Ec. (17) es exacta con \(p=p(n)\)',
     r'Ecs. (A2)–(A4): transformación de Liouville y potencial \(\mathcal{U}\)',
     r'Ec. (A6): ecuación de Milne para la amplitud \(Q\)',
     r'Ecs. (A7) y (A18): \(\mathcal{L}_\pm\) y su inversión para \((Q, Q^\prime)\)',
     r'Ec. (A10): \(\mathcal{L}_\mathrm{ext}\) a partir de \(\mathcal{L}_+\)',
     r'Ec. (A13): \(\mathcal{U}\) hasta \(\mathcal{O}(\tilde r^{-2})\)',
     r'Ec. (A14): reducción a la ecuación de Whittaker (\(z,\kappa,\mu\))',
     r'Ecs. (21)–(22): forma de Schrödinger en coordenada tortuga y \(\Delta V\)',
     r'Ec. (23): niveles hidrogénicos (\(n=l+1,\ l+2\))',
     r'Ec. (B4): flujo \(\sqrt{-g}\,T^r{}_t\) del axión',
     r'Ec. (B6): luminosidad promediada en una oscilación',
     r'Ec. (B7b): energía del modo en el axión',
     r'Ec. (B7a): energía del modo en el fluido, desde \(T^t{}_t\)'],
    [r'\(\lambda_a = 19.7\) km para \(m_a = 10^{-11}\) eV (paper: \(\simeq 20\) km)',
     r'\(\nu_a = 2.42\) kHz para \(m_a = 10^{-11}\) eV (paper: \(\simeq 2.4\) kHz)',
     r'Tabla I: desplazamientos de 18.1 %, 14.5 % y 19.0 % (texto: 18, 15 y 19 %)',
     r'Las cifras del texto coinciden con las Tablas I y II; el acoplo compensa el 84 % del desplazamiento de \(f\)',
     r'Tabla II: \(\mathrm{Im}\,\nu \neq 0\) exactamente cuando \(\mathrm{Re}\,\nu > \nu_a\) (36 entradas)',
     r'\(\tau_a\) mínimo \(= 12\) ms (texto: \(\sim 10\) ms)',
     r'Tabla III: el asterisco y la Fig. 4 coinciden con la tabla',
     r'Tabla IV: \(\alpha = GMm_a\) reproduce la columna (dif. máx. 0.0005)',
     r'Tabla IV: \(\nu_n\) coherente con \(-\tilde k_n^2 n^2/\alpha^2\)',
     r'Tabla IV recalculada desde el exterior: hasta 10 kHz, dif. máx. 5 Hz',
     r'Tabla IV, 11–15 kHz: valores del paper dentro de la banda que deja el interior',
     r'Cruces evitados: niveles \(n=3\) (12.19 kHz) y \(n=4\) (12.06 kHz) en 11.92–11.94 kHz, frente a \(p_2\) = 11.93 kHz',
     r'Figura generada'],
    # ---- Paso 7: más allá de Cowling (trabajo nuevo)
    [r'Parte sin traza de las ecuaciones angulares: \(H_2 = H_0\), sin contribución del axión',
     r'Con la métrica congelada se recuperan exactamente las ecuaciones de Cowling del paso 5',
     r'Bianchi: \(G_{tt}\) linealizada es redundante',
     r'Bianchi: la traza angular linealizada es redundante',
     r'La ligadura algebraica de \(H_0\) se conserva al integrar',
     r"Sin axión: \(H_1'\), \(K'\), \(W'\), \(H_0\) y \(V\) de Lindblom-Detweiler (versión de Zhao y Lattimer)",
     (r"\(X'\) de Zhao y Lattimer: errata de signo en el término \(2\omega^2 W\)", 'obs'),
     (r'Burgio et al.: cinco de las seis ecuaciones tienen erratas de transcripción', 'obs'),
     r"Kunjipurayil et al. (2022, coautor T. Zhao): las seis coinciden, con \(-2\omega^2\) en \(X'\)",
     (r'Zhang et al. (2026): escritas ya con las cinco correcciones a Burgio et al.; errata propia, \(\omega^3\) en \(V\)', 'obs'),
     r"Control negativo: \(16\pi\to8\pi\) en \(H_1'\) se detecta",
     r'Exterior sin axión: la función de Zerilli cumple su ecuación',
     r'Exterior con axión: \(\delta\theta\) se desacopla de la métrica y cumple la Ec. (13) del Paper II',
     r'Centro: \(H_0(0)\), \(H_1(0)\) y \(V(0)\) como en LD, sin corrección axiónica',
     r'Centro: \(X(0)\) con la corrección axiónica',
     r"Centro: \(\delta\theta''(0)\) = Ec. (16c) \(-\,3K(0)\theta_0''(0)/(2l+3)\)",
     r"Sistema extendido: \(H_1'\)",
     r"Sistema extendido: \(K'\)",
     r"Sistema extendido: \(W'\), sin cambios",
     r"Sistema extendido: \(X'\)",
     r'Sistema extendido: relación algebraica de \(H_0\)',
     r'Sistema extendido: relación algebraica de \(V\)',
     r'Sistema extendido: Klein-Gordon, Ec. (8d) + términos métricos',
     r'Todos los términos nuevos se anulan sin axión',
     r"Control negativo: quitar el flujo axiónico de \(K'\) se detecta"],
    # ---- Paso 8: números y figuras de la tercera parte
    [r'Fondo: el modelo de 1.4 \(M_\odot\) de Lindblom, Mendell e Ipser sale con \(R = 14.45\) km',
     r'Cowling (Ecs. 8a–8c): los 24 modos \(f, p_1, p_2\) de Yoshida y Kojima, 23 dentro de su última cifra y uno a 0.0057',
     r'Métrica dinámica sin axión: Re \(\omega\) de los modos \(f\) con \(l=2\)–5, igual a las 5 cifras de Lindblom, Mendell e Ipser',
     r'Im \(\omega\) (emisión de ondas gravitacionales) dentro del 2 % de Lindblom, Mendell e Ipser',
     r'Métrica dinámica: los 24 modos «exactos» de Yoshida y Kojima, dentro del 0.5 %',
     r'Convergencia numérica: Re \(\omega\) estable a \(10^{-5}\)',
     (r"Con el \(X'\) impreso por Zhao y Lattimer, el modo \(f\) sale un 31 % más bajo", 'obs'),
     r'Fondo con axión, \(\nu_a = 5\) kHz: converge, \(\theta_0 \approx \pi\) dentro y \(\to 0\) fuera',
     r'Fondo con axión, \(\nu_a = 15\) kHz: converge, \(\theta_0 \approx \pi\) dentro y \(\to 0\) fuera'],
    # ---- Paso 9: notas con axión en Cowling sobre el fondo BSk26 (Tablas I, II y IV del Paper II)
    [r'Fondo BSk26 sin axión, \(p_c = 100\) MeV/fm³: \(M = 1.536\,M_\odot\), \(R = 11.72\) km (Paper I: 1.53 y 11.7)',
     r'Cowling sin axión: \(f, p_1, p_2\) = 2.429, 7.171, 10.194 kHz frente a 2.43, 7.15, 10.18 (máx. 0.29 %)',
     r'Umbral del condensado: no hay a 2.5 kHz y sí a 2.8 kHz, como en el Paper I',
     r'Tabla IV: masas de las estrellas con axión de 8 a 15 kHz, dentro del {{V:c9_t4}}',
     r'Tabla I (solo fondo): los 24 desplazamientos \(\nu/\nu_0-1\), dentro de {{V:c9_t1}} puntos',
     r'Tabla II (acoplado), parte real: los 36 desplazamientos dentro de {{V:c9_t2}} puntos; el modo \(f\), de {{V:c9_t2f}}',
     ('{{V:c9_lbl6}}', '{{K:c9_lbl6}}'),
     r'Tabla II, amortiguamiento: los {{V:c9_nratio}} valores \(\geq 1\) Hz dentro de un factor 2; \(\tau_a\) mínimo {{V:c9_tau}} (paper: 12 ms)'],
    # ---- Paso 10: con axión y métrica dinámica (etapa B)
    [r'Conversión del sistema extendido a números: sin axión, igual que la del paso 8 (\(10^{-15}\))',
     r'Sustituciones del axión: las 8 cantidades del Paper II (\(n\,\partial m^*/\partial\theta\), \(m_\mathrm{eff}^2\), \(\varepsilon(n,\theta)\), \(c_s^2\)…)',
     r'Todo el programa, sin axión, con el politropo de LMI: Re \(\omega\) = 1.20162 (LMI: 1.2016), Im \(\omega\) = 3.944\(\times10^{-4}\) (LMI: 3.95\(\times10^{-4}\))',
     r'Convergencia del modo \(f\) a 6 kHz: Re \(\nu\) estable a {{V:c10_conv_re}}, Im \(\nu\) al {{V:c10_conv_im}}',
     r'Sin axión, la métrica dinámica baja las notas de Cowling: \(f, p_1, p_2\) = 1.961, 6.594, 9.872 kHz (−19, −8, −3 %); \(\tau_f\) = 0.174 s',
     r'Las 36 notas acopladas, encontradas y amortiguadas',
     r'Umbral: los modos \(p\) radian axiones solo con \(\mathrm{Re}\,\nu > \nu_a\)',
     (r'Con la métrica dinámica, el axión sube \(f\) hasta un {{V:c10_fmax}} (Cowling: {{V:c9_fmax}}) y acorta su amortiguamiento de {{V:c10_tau0}} a {{V:c10_tau15}}', 'obs'),
     (r'Amortiguamiento máximo: \(p_2\) a {{V:c10_p2nu}}, {{V:c10_p2}} (Cowling: {{V:c9_p2}}); la suma \(1/\tau_{GW}+1/\tau_a\) solo acierta en orden de magnitud', 'obs')],
    # ---- Pasos 14 y 15: modos radiales (informe aparte); se muestran los nombres que genera el propio paso
    None,
    None,
    None,
    None,
    None,
    None,
    None,
]

RADIAL_URL = 'https://claude.ai/artifact/9XoYYWX9wvNUKoav2HmQbB'
BASE = Path(__file__).parent / 'baseline_v1'     # datos de la versión publicada antes de corregir el lapso
BASE_DATE = '2026-11-05'


# ---------------------------------------------------------------- cifras del texto, generadas desde los datos
def _rx(text, pat):
    m = re.search(pat, text)
    assert m, (pat, text)
    return m.group(1)


def values(checks, axcow, axgr):
    # Cada cifra del texto de las partes 4 y 5 sale de checks.json o de las figuras; se calcula igual para la
    # versión anterior (baseline_v1) y la actual, y las que difieren se marcan en el HTML.
    by = {s['module']: [c['name'] for c in s['checks']] for s in checks['steps']}
    c9, c10 = by['step9_axion_cowling'], by['step10_axion_gr']
    nus = axcow['nus']
    assert nus == axgr['nus'] and nus[-1] == 15.0
    cf, cp1, cp2 = (m for m in axcow['modes'])
    gf, gp1, gp2 = (m for m in axgr['modes'])
    at = lambda arr, nu: arr[nus.index(nu)]
    V = {}
    V['c9_t4'] = _rx(c9[3], r'máx\. ([\d.]+ %)')
    V['c9_t1'] = _rx(c9[4], r'máx\. ([\d.]+)\)')
    V['c9_t2'] = _rx(c9[5], r'máx\. ([\d.]+);')
    V['c9_t2f'] = _rx(c9[5], r'modo f, ([\d.]+)\)')
    V['c9_tmax'] = max(V['c9_t1'], V['c9_t2'], key=float)
    V['c9_nratio'] = _rx(c9[7], r'los (\d+) \|Im')
    V['c9_ratio'] = _rx(c9[7], r'cociente ([\d.]+-[\d.]+)').replace('-', '–')
    V['c9_tau'] = _rx(c9[7], r'tau_a mínimo = ([\d.]+ ms)')
    V['c9_fmax'] = f"{max(cf['cp']):.1f} %"
    V['c9_f15'] = f"{cf['cp'][-1]:.1f} %"
    # p1 a 8 kHz en Cowling: ¿radia? (Re nu = nu_0 (1 + desplazamiento))
    re8 = axgr['cow0'][1] * (1 + at(cp1['cp'], 8.0) / 100)
    gap8, rad8 = 100 * (re8 / 8.0 - 1), at(cp1['im'], 8.0) > 1e-3
    exc = 'sin excepciones' not in c9[6]
    assert exc == rad8, (c9[6], rad8)
    if rad8:
        V['c9_edge'] = (rf'Los modos que radian son los mismos, salvo \(p_1\) a 8 kHz, cuya parte real cae a un '
                        rf'{abs(gap8):.1f} % del umbral: en nuestro cálculo queda justo por encima.')
        V['c9_lbl6'] = r'Tabla II: los mismos modos amortiguados, salvo \(p_1\) a 8 kHz, que cae justo en el umbral'
    else:
        V['c9_edge'] = (rf'Los modos que radian son los mismos, sin excepción; el caso más ajustado es \(p_1\) a 8 kHz, '
                        rf'cuya parte real queda un {abs(gap8):.2f} % por debajo del umbral.')
        V['c9_lbl6'] = r'Tabla II: los mismos modos amortiguados que en la secuela, sin excepciones'
    V['K:c9_lbl6'] = 'obs' if rad8 else ''
    i5 = max(range(len(nus)), key=lambda i: gf['cp'][i])
    V['c10_fmax'] = f"{gf['cp'][i5]:.1f} %"
    assert V['c10_fmax'] == _rx(c10[7], r'hasta un ([\d.]+ %)'), c10[7]
    V['c10_fmax_nu'] = f'{nus[i5]:.0f} kHz'
    V['c10_f15'] = f"{gf['cp'][-1]:.1f} %"
    V['c10_f6'] = f"{at(gf['re'], 6.0):.3f} kHz"
    m = re.search(r'Re nu en (\d)e-0?(\d+) e Im nu en ([\d.]+ %)', c10[3])
    assert m, c10[3]
    V['c10_conv_re'] = rf'\({m.group(1)}\times10^{{-{m.group(2)}}}\)'
    V['c10_conv_im'] = m.group(3)
    V['c10_tau0'] = _rx(c10[4], r'gravitacionales en ([\d.]+ s)')
    V['c10_tau15'] = _rx(c10[7], r'pasa de [\d.]+ a ([\d.]+ s)')
    V['c10_pdiff'] = f"{max(abs(g - c) for G, C in ((gp1, cp1), (gp2, cp2)) for g, c in zip(G['cp'], C['cp'])):.1f}"
    gre8 = at(gp1['re'], 8.0)
    if rad8:
        V['c10_edge'] = r'\(p_1\) deja de radiar a 8 kHz, donde en nuestro Cowling (paso 9) aún radia un poco.'
    else:
        V['c10_edge'] = (rf'A 8 kHz, \(p_1\) ({gre8:.2f} kHz) queda un {100 * (1 - gre8 / 8):.0f} % por debajo de \(\nu_a\) '
                         rf'y solo radia ondas gravitacionales; en Cowling (paso 9) tampoco radia, aunque allí queda a un '
                         rf'{abs(gap8):.2f} % del umbral.')
    j = max(range(len(nus)), key=lambda i: gp2['im'][i])
    V['c10_p2nu'], V['c10_p2'] = f'{nus[j]:.0f} kHz', f"{gp2['im'][j]:.1f} Hz"
    assert f"p2 a {nus[j]:.0f} kHz, {gp2['im'][j]:.1f} Hz" in c10[8], c10[8]
    V['c10_p2tau'] = f"{1e3 / (2 * math.pi * gp2['im'][j]):.1f} ms"
    V['c9_p2'] = f"{cp2['im'][j]:.1f} Hz"
    V['c10_p2_8'], V['c9_p2_8'] = f"{at(gp2['im'], 8.0):.2f} Hz", f"{at(cp2['im'], 8.0):.1f} Hz"
    return V


DESCR = {
    'c9_t4': 'Tabla IV: diferencia máxima en las masas (Cowling)',
    'c9_t1': 'Tabla I: diferencia máxima en los desplazamientos (puntos)',
    'c9_t2': 'Tabla II: diferencia máxima en los desplazamientos (puntos)',
    'c9_t2f': 'Tabla II, modo f: diferencia máxima (puntos)',
    'c9_tmax': 'Tablas I y II: diferencia máxima (puntos)',
    'c9_nratio': 'Tabla II: amortiguamientos ≥ 1 Hz comparados',
    'c9_ratio': 'Tabla II: cociente de amortiguamientos (este cálculo / secuela)',
    'c9_tau': 'Cowling: τ_a mínimo',
    'c9_fmax': 'Cowling: subida máxima del modo f',
    'c9_f15': 'Cowling: subida del modo f a 15 kHz',
    'c9_edge': 'Cowling: p₁ a 8 kHz frente al umbral',
    'c9_lbl6': 'Comprobación: modos amortiguados frente a la secuela',
    'c10_fmax': 'Métrica dinámica: subida máxima del modo f',
    'c10_fmax_nu': 'Métrica dinámica: ν_a de la subida máxima de f',
    'c10_f15': 'Métrica dinámica: subida del modo f a 15 kHz',
    'c10_f6': 'Métrica dinámica: modo f a 6 kHz',
    'c10_conv_re': 'Métrica dinámica: estabilidad de Re ν del modo f (6 kHz)',
    'c10_conv_im': 'Métrica dinámica: estabilidad de Im ν del modo f (6 kHz)',
    'c10_tau0': 'Métrica dinámica: τ_GW del modo f sin axión',
    'c10_tau15': 'Métrica dinámica: τ_GW del modo f a 15 kHz',
    'c10_pdiff': 'Modos p: diferencia máxima con Cowling (puntos)',
    'c10_edge': 'Métrica dinámica: p₁ a 8 kHz frente al umbral',
    'c10_p2nu': 'Métrica dinámica: ν_a del amortiguamiento máximo de p₂',
    'c10_p2': 'Métrica dinámica: amortiguamiento máximo de p₂',
    'c10_p2tau': 'Métrica dinámica: τ_a de p₂ en ese máximo',
    'c9_p2': 'Cowling: amortiguamiento de p₂ a esa ν_a',
    'c10_p2_8': 'Métrica dinámica: amortiguamiento de p₂ a 8 kHz',
    'c9_p2_8': 'Cowling: amortiguamiento de p₂ a 8 kHz',
}


def _mark(new, old):
    if old is None or new == old:
        return new
    return f'<mark class="nuevo" title="Antes: {html.escape(old)}">{new}<span class="antes">{old}</span></mark>'


def fill(text, V, V0, plain=False):
    def sub(m):
        key = m.group(1)
        return V[key] if plain else _mark(V[key], V0.get(key))
    return re.sub(r'\{\{V:(\w+)\}\}', sub, text)


def _decimals(text):
    # cifras decimales de un texto, sin LaTeX en línea ni números de ecuación («Ec. (4.13)», «(4.1)»)
    text = re.sub(r'\\\(.*?\\\)', ' ', text)
    text = re.sub(r'Ecs?\.\s*\(?[\d.]+[a-z]?|\(\d+\.\d+[a-z]?\)', ' ', text)
    return set(re.findall(r'(?<![\w.])\d+\.\d+', text))


def cambios(V, V0, checks, base_checks, figs_changed):
    keys = [k for k in DESCR if V0.get(k) != V[k] and k not in ('c9_edge', 'c10_edge', 'c9_lbl6', 'c9_tmax')]
    rows_ = ''.join(f'<tr><td>{html.escape(DESCR[k])}</td><td class="num">{V0[k]}</td><td class="num"><b>{V[k]}</b></td></tr>'
                    for k in keys)
    base = {s['module']: s['checks'] for s in base_checks['steps']}
    nch = sum(1 for s in checks['steps'] if s['module'] in base
              for c, c0 in zip(s['checks'], base[s['module']]) if c['name'] != c0['name'])
    flips = [c['name'] for s in checks['steps'] if s['module'] in base
             for c, c0 in zip(s['checks'], base[s['module']]) if c['ok'] != c0['ok']]
    nflip_txt = 'ninguna cambia de resultado' if not flips else f'{len(flips)} cambian de resultado'
    new_steps = [s['step'] for s in checks['steps'] if s['module'] not in base]
    better = float(V['c9_tmax']) < float(V0['c9_tmax'])
    return f'''<section class="stack cambios" id="cambios">
    <p class="eyebrow">Versión 2 · {checks['env']['date']} · sustituye a la del {BASE_DATE}</p>
    <h2>Qué cambia en esta versión</h2>
    <p><b>Corrección en el código compartido de la estrella con axión.</b> El lapso se normalizaba con \\(\\Phi(\\infty)=0\\) imponiendo en la superficie \\(e^{{2\\Phi(R)}}=1-2M(R)/R\\), sin contar la masa de la cola del condensado que queda fuera de la estrella. Ahora la condición se impone al final de la integración, donde la cola ya ha decaído. Afecta a todo lo calculado con axión (cuarta y quinta parte); sin axión no cambia nada. Lo destapó una revisión externa del informe de modos radiales, en el que el resultado dependía del radio de empalme.</p>
    <p>Las cifras que cambian aparecen así: <mark class="nuevo">valor nuevo<span class="antes">valor anterior</span></mark>. Las gráficas de la cuarta y quinta parte están recalculadas{'' if figs_changed else ' (sin cambios visibles)'}. {nch} comprobaciones cambian de cifras y llevan la etiqueta «actualizada»; {nflip_txt}. El acuerdo con las Tablas I y II de la secuela {'mejora' if better else 'no mejora'}: la diferencia máxima pasa de {V0['c9_tmax']} a {V['c9_tmax']} puntos.{(' ' + _mark(V['c9_edge'], V0['c9_edge'])) if V['c9_edge'] != V0['c9_edge'] else ''}</p>
    <p>Se añaden {len(new_steps)} pasos ({', '.join(html.escape(x) for x in new_steps)}), con <a href="{RADIAL_URL}">informe propio</a>, y se retiran los enlaces a las tres «líneas abiertas», cuyo informe no resistió la auditoría.</p>
    <div class="tw"><table><thead><tr><th>Cifra</th><th>Antes</th><th>Ahora</th></tr></thead><tbody>{rows_}</tbody></table></div>
  </section>'''



def rows(V, V0):
    data = json.loads((RES / 'checks.json').read_text())
    base = {st['module']: st['checks'] for st in json.loads((BASE / 'checks.json').read_text())['steps']}
    assert len(data['steps']) == len(PRETTY), (len(data['steps']), len(PRETTY))
    out, n_ok, n_tot, per_step = [], 0, 0, []
    for step, pretty in zip(data['steps'], PRETTY):
        checks = step['checks']
        old = base.get(step['module'])
        if pretty is None:                                  # pasos sin etiqueta a mano: el nombre generado
            pretty = [html.escape(c['name']) for c in checks]
        assert len(checks) == len(pretty), (step['step'], len(checks), len(pretty))
        items = []
        for i, (c, p) in enumerate(zip(checks, pretty)):
            label, kind = (p if isinstance(p, tuple) else (p, ''))
            kind = re.sub(r'\{\{K:(\w+)\}\}', lambda m: V['K:' + m.group(1)], kind)
            # guarda: toda cifra decimal tecleada en la etiqueta debe estar en el nombre que genera el paso
            missing = _decimals(fill(label, V, V0, plain=True)) - _decimals(c['name'])
            assert not missing or '{{V:' in label, (step['module'], i, sorted(missing), c['name'])
            label = fill(label, V, V0)
            n_tot += 1
            n_ok += c['ok']
            if not c['ok']:
                chip = '<span class="chip chip-bad">no cuadra</span>'
            elif kind == 'obs':
                chip = '<span class="chip chip-obs">observación</span>'
            else:
                chip = '<span class="chip chip-ok">coincide</span>'
            if old is not None and c['name'] != old[i]['name']:
                chip = f'<span class="chips">{chip}<span class="chip chip-new">actualizada</span></span>'
            items.append(f'<li title="{html.escape(c["name"])}">{chip}<span class="ck">{label}</span></li>')
        tag = '' if old is not None else ' <span class="chip chip-new">nuevo en esta versión</span>'
        out.append(f'<section class="steplist"><h3>{html.escape(step["step"])}{tag}</h3><ul class="checks">'
                   + ''.join(items) + '</ul></section>')
        per_step.append((sum(c['ok'] for c in checks), len(checks)))
    return '\n'.join(out), n_ok, n_tot, data, per_step


def repro(data, per_step):
    steps = ''.join(
        f'<tr><td>{i}</td><td>{html.escape(s["step"])}</td><td class="mod"><code>{s["module"]}</code></td>'
        f'<td class="num">{ok}/{tot}</td><td class="num">{s["seconds"]:.0f} s</td></tr>'
        for i, (s, (ok, tot)) in enumerate(zip(data['steps'], per_step), 1))
    env = data['env']
    envrows = [('Python', env['python']), ('SymPy', env['sympy']), ('mpmath', env['mpmath']),
               ('NumPy', env['numpy']), ('SciPy', env['scipy']), ('Matplotlib', env['matplotlib']),
               ('Sistema', env['platform']), ('CPU', env['cpu']),
               ('Última ejecución', env['date'])]
    return steps, ''.join(f'<dt>{k}</dt><dd>{html.escape(v)}</dd>' for k, v in envrows)


def main():
    tpl = (Path(__file__).parent / 'template.html').read_text()
    checks = json.loads((RES / 'checks.json').read_text())
    load = lambda d, n: json.loads((d / f'fig_{n}.json').read_text())
    V = values(checks, load(RES, 'axcow'), load(RES, 'axgr'))
    V0 = values(json.loads((BASE / 'checks.json').read_text()), load(BASE, 'axcow'), load(BASE, 'axgr'))
    figs_changed = any(load(RES, n) != load(BASE, n) for n in ('axcow', 'axgr'))
    table, n_ok, n_tot, data, per_step = rows(V, V0)
    steps, env = repro(data, per_step)
    fails = '\n'.join('  FAIL: ' + html.escape(c['name']) for s in data['steps'] for c in s['checks'] if not c['ok'])
    ok1, tot1 = (sum(x[i] for x in per_step[0:4]) for i in (0, 1))
    ok2, tot2 = (sum(x[i] for x in per_step[4:6]) for i in (0, 1))
    ok3, tot3 = per_step[6]
    ok4, tot4 = per_step[7]
    ok5, tot5 = per_step[8]
    ok6, tot6 = per_step[9]
    ok7, tot7 = (sum(x[i] for x in per_step[10:17]) for i in (0, 1))
    figdata = json.dumps({k: json.loads((RES / f'fig_{f}.json').read_text())
                          for k, f in [('toy', 'toy'), ('levels', 'levels'), ('modes', 'modes'), ('meff', 'meff'),
                                       ('background', 'background'), ('axcow', 'axcow'), ('axgr', 'axgr')]}, separators=(',', ':')).replace('</', '<\\/')
    charts_js = (Path(__file__).parent / 'charts.js').read_text()
    anim_js = (Path(__file__).parent / 'anim.js').read_text()
    page = (tpl.replace('{{CHECKS}}', table).replace('{{N_OK}}', str(n_ok)).replace('{{N_TOT}}', str(n_tot))
            .replace('{{N_OK1}}', str(ok1)).replace('{{N_TOT1}}', str(tot1))
            .replace('{{N_OK2}}', str(ok2)).replace('{{N_TOT2}}', str(tot2))
            .replace('{{N_OK3}}', str(ok3)).replace('{{N_TOT3}}', str(tot3))
            .replace('{{N_OK4}}', str(ok4)).replace('{{N_TOT4}}', str(tot4))
            .replace('{{N_OK5}}', str(ok5)).replace('{{N_TOT5}}', str(tot5))
            .replace('{{N_OK6}}', str(ok6)).replace('{{N_TOT6}}', str(tot6))
            .replace('{{N_OK7}}', str(ok7)).replace('{{N_TOT7}}', str(tot7))
            .replace('{{CAMBIOS}}', cambios(V, V0, checks, json.loads((BASE / 'checks.json').read_text()), figs_changed))
            .replace('{{RADIAL_URL}}', RADIAL_URL)
            .replace('{{FIGDATA}}', figdata).replace('{{CHARTS_JS}}', charts_js).replace('{{ANIM_JS}}', anim_js)
            .replace('{{SECS}}', str(data['seconds'])).replace('{{MINS}}', str(round(data['seconds'] / 60))).replace('{{STEPS}}', steps).replace('{{ENV}}', env)
            .replace('{{FAILS}}', fails))
    page = fill(page, V, V0)
    assert '{{' not in page, re.findall(r'\{\{[^}]*\}\}', page)[:5]
    out = RES / 'informe_axiones.html'
    out.write_text(page)
    guia_out = RES / 'guia_departamento.html'
    guia_out.write_text((Path(__file__).parent / 'guia_departamento.html').read_text())
    print(out, f'{len(page) / 1e3:.0f} kB', f'{n_ok}/{n_tot}')
    print(guia_out, f'{guia_out.stat().st_size / 1e3:.0f} kB')
    for extra in ():                       # las tres «líneas» retiradas ya no se copian (LAB_NOTEBOOK §3)
        src = Path(__file__).parent / extra
        if src.exists():
            dst = RES / extra
            dst.write_text(src.read_text())
            print(dst, f'{dst.stat().st_size / 1e3:.0f} kB')


if __name__ == '__main__':
    main()
