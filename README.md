# Workshop física + IA — caso 2: cálculo simbólico

Reproducción, desde la acción, de las ecuaciones de dos papers:

> A. Gómez-Bañón, P. Pnigouras, J. A. Pons, *Axion condensates in neutron stars and radial
> oscillation modes*, [arXiv:2510.16245](https://arxiv.org/abs/2510.16245) (Paper I).
>
> A. Gómez-Bañón, P. Pnigouras, J. A. Pons, *Nonradial oscillation modes of axion condensed neutron
> stars in the Cowling approximation*, [arXiv:2610.37308](https://arxiv.org/abs/2610.37308) (Paper II).

El flujo es el típico de Mathematica hecho a mano, pero automatizado con SymPy y con
**verificación en cada paso** contra las ecuaciones de los papers, transcritas desde su LaTeX en
`axion_ns/paper.py` y `axion_ns/nonradial.py`. La numeración de ecuaciones es la de los PDF publicados.
El Paper I numera por secciones, (2.6), (3.3)…; el Paper II, de forma secuencial, (8a), (13), (A6)…

```
acción  S = ∫√-g [ -R/16πG + ½(∂a)² - V(a) - ε(n,a) ],   ε = ε₀(n) + n (m*(a) - m_N)
  │  Einstein (G_μν y T_μν por variación métrica) + Klein-Gordon (Euler-Lagrange)
  │  + conservación de T y del número bariónico
  ▼
Paper I
paso 1  fondo estático  →  TOV + KG, Ec. (2.6a-d)
paso 2  perturbaciones radiales lineales  →  Ec. (3.3a-e) y exterior Ec. (3.7)
paso 3  modelo simplificado (Sec. IV)  →  Ec. (4.2)-(4.13): matriz, dispersión, ec. trascendente
paso 4  cifras numéricas del texto + reproducción de la Fig. 3 (QNMs del modelo simplificado)
Paper II
paso 5  no radiales en Cowling  →  Ec. (8a-d), (10), (13), (15)-(17), (21)-(23), apéndices A y B
paso 6  tablas I-IV + cálculo independiente de los niveles ligados (Tabla IV, Fig. 5)
Trabajo nuevo
paso 7  más allá de Cowling: métrica perturbada → Lindblom-Detweiler + axión, centro y exterior
paso 8  números: modos f/p de un politropo frente a la literatura; fondo con axión; figuras
paso 9  con axión en Cowling: estrella BSk26 con condensado; Tablas I, II y IV del Paper II
paso 10 con axión y métrica dinámica (sin publicar): f, p1, p2 para 12 masas del axión
```

## Ejecutar

```bash
uv sync --frozen                              # entorno con las versiones de uv.lock (Python en .python-version)
uv run python -m axion_ns.run_all             # todo (~75 min); esperado: 128/131, con 3 FAIL que son las cifras del Paper I
uv run python -m axion_ns.step5_nonradial     # un paso suelto
uv run python report/build.py                 # regenera results/informe_axiones.html y results/axion_ns_codigo.zip
```

Las figuras salen en `results/toy_qnm_landscape.png` y `results/paper2_atom_levels.png`.

## Cómo se verifica

- Las ecuaciones de los papers se comparan **módulo las ecuaciones de fondo**: se despeja una
  variable de la ecuación derivada, se sustituye en la del paper, se eliminan derivadas del fondo
  con TOV+KG y se simplifica hasta 0.
- **Controles negativos**: se plantan errores (4π→8π, un signo, un término omitido, 2l por 2(l+1))
  y el comparador tiene que dar FAIL.
- **Consistencia interna**: identidad ∇T(axión) = (□a+V′)∂a, redundancia de G_tt linealizada
  (Bianchi), límite sin axión → TOV estándar, exterior → KG en Schwarzschild.
- Paper II: la dependencia angular se elimina con la ecuación de Legendre (modo m = 0; las
  ecuaciones no dependen de m), y las condiciones en el origen salen de un desarrollo en serie con
  un fondo regular genérico.

## Resultado Paper I (40 comprobaciones, 37 OK)

**Todo el desarrollo analítico del paper es correcto**: Ec. (2.6), (3.3), (3.7) y toda la cadena
del modelo simplificado (4.2)→(4.13), incluida la ecuación trascendente obtenida casando interior y
exterior con δp(x₀)=0. Los modos de la Fig. 3 se reproducen.

Observaciones para comentar con los autores:

1. **Punto de linealización del modelo simplificado (Sec. IV).** Las aproximaciones (4.1),
   ∂m*/∂a ≈ −βσ_N/4f_a y V″+n m*″ ≈ m_a²(1−α), corresponden a θ = a/f_a ≈ 0 (y la primera, de hecho,
   a sin θ/g ≈ 1). Pero el texto toma θ = π en el interior cuando α>1. En θ = π:
   - ∂m*/∂a = −βσ_N sin θ /(4 f_a g) = **0**: no hay acoplo lineal fluido-axión;
   - m_eff² = m_a²(α−1)/√(1−β) ≈ 26 m_a² para α=10, **positiva**, frente a 1−α = −9 (taquiónica).

   Linealizando de forma consistente en θ = π, los modos fluidos quedan sin amortiguar y el espectro
   axiónico cambia (panel derecho de la figura). En la estrella realista el acoplo vendría de las
   zonas donde a₀ no está en un mínimo, como la capa de transición. Las secciones V y siguientes usan
   m*(a) y V(a) completos, así que no les afecta. El Paper II describe justamente ese fondo
   (θ∼π dentro, acoplo solo en la capa de transición; pie de la Fig. 5 y Sec. III E).
2. **m_a R ≃ 0.12** para 3 kHz y 12 km: sale de νR/c. En unidades naturales m_a R = 2πνR/c ≈ **0.75**
   (coherente con 3 kHz ↔ 12 peV, E = hν). La longitud de Compton es ≈1.3 R, no "varios radios".
3. **m_a² f_a² ≃ ε·13 MeV/fm³**: con m_a² f_a² = εβ m_π² f_π²/4 sale ε·(4.4–4.8) MeV/fm³. El
   Paper II usa ε m_π² f_π² ≃ 2 MeV/fm³ para ε = 0.1 (pie de la Fig. 7), que es coherente con
   esto. Además, **λ² ≃ ε·0.02** no se sigue ni del 13 (13/150 = 0.087) ni del valor recalculado
   (≈0.03).
4. Menores: la Ec. (2.8), p(n,a) ≈ p(n,0), es **exacta** dada la Ec. (2.7). Y en la forma de la
   Ec. (4.13) con tangentes cada modo fluido es un cero pegado a un polo de tan(k₊x₀). Para buscar
   raíces conviene la forma entera f·cos(k₊x₀)cos(k₋x₀), que parece ser la representada en la Fig. 3.

## Resultado Paper II (40 comprobaciones, 40 OK)

- Las ecuaciones de perturbación (8a-d) salen de la acción, igual que el exterior (13), el ajuste
  (15), las condiciones en el origen (10a,b), (16a-c), (17), toda la maquinaria del apéndice A
  (Liouville, Milne, Whittaker), la forma de Schrödinger (21)-(22), los niveles (23) y los flujos y
  energías del apéndice B.
- Las tablas I–IV son coherentes entre sí y con el texto: el umbral Im ν ≠ 0 ⇔ Re ν > ν_a se
  cumple en las 36 entradas de la Tabla II, y τ_a mínimo = 12 ms.
- **Tabla IV reproducida de forma independiente.** Se resuelve solo la Ec. (13) en Schwarzschild
  con la condición de Whittaker en el infinito; el interior entra únicamente como condición en R,
  acotada entre "δθ plano" y Dirichlet, con R = 10.5–12 km. Hasta ν_a = 10 kHz, los niveles n = 3, 4
  coinciden con la tabla a ≤5 Hz con cualquier interior. Por encima, el interior pesa (banda de
  56 Hz a 15 kHz) y los valores del paper caen dentro. Los cruces evitados con p₂ (12.19 y
  12.06 kHz) salen en 11.92–11.94 kHz, frente a los 11.93 kHz del texto.

## Paso 7: más allá de Cowling (trabajo nuevo)

Es el primer "natural next step" de la Sec. IV del Paper II. Se perturba también la métrica (gauge de
Regge-Wheeler, variables de Lindblom y Detweiler) y se deriva el sistema completo desde la misma
acción. No hay ecuaciones publicadas con las que comparar, así que se valida así:

- la parte sin traza angular da H2 = H0; G_tt y la traza angular son redundantes (Bianchi) y la
  ligadura algebraica de H0 se conserva;
- con la métrica congelada salen exactamente las (8a)-(8d) del paso 5;
- sin axión salen las ecuaciones de LD, contrastadas con cuatro transcripciones: Zhao y Lattimer
  (arXiv:2204.03037) tienen un signo cambiado en X', ya corregido en Kunjipurayil et al.
  (arXiv:2205.02081); Burgio et al. (arXiv:1106.2736) tienen erratas en cinco de seis ecuaciones, ya
  corregidas en Zhang et al. (arXiv:2601.21911), que a su vez ponen w^3 en V;
- en el exterior, Zerilli para la métrica y la Ec. (13) del Paper II para el axión, desacoplados;
- en el centro, las condiciones de LD más las correcciones axiónicas en X(0) y δθ''(0).

El sistema extendido, con los términos nuevos separados, está en `axion_ns/fullgr.py` (`NEW`, `EXT`)
y en la tercera parte de la página. Con axión no está integrado numéricamente.

## Paso 8: comprobación numérica y figuras

`axion_ns/fmodes.py` integra las ecuaciones derivadas sin axión (fondo TOV, condiciones de LD en el
centro, X(R) = 0, Zerilli fuera; el modo es el mínimo de la amplitud saliente a frecuencia real):

- el sistema extendido reproduce las 5 cifras de Re w de los modos f l = 2-5 del politropo de
  Lindblom, Mendell e Ipser (1997, gr-qc/9704046, Tabla I) e Im w dentro del 2 %;
- las (8a)-(8c) de la secuela reproducen los 24 modos de Cowling de Yoshida y Kojima (1997,
  gr-qc/9705081): 23 dentro de su última cifra y uno a 0.0057; sus 24 valores exactos, dentro del 0.5 %;
- con el signo de X' impreso por Zhao y Lattimer, el modo f sale un 31 % más bajo;
- el fondo con axión (politropo + Klein-Gordon, sin retroacción) muestra dónde actúan los términos
  nuevos.

Ojo: las dos referencias usan p = K eps^2 con eps la densidad de energía.

## Paso 9: con axión, en Cowling (etapa A)

`axion_ns/bsk26.py` implementa BSk26 con los ajustes analíticos de Potekhin et al. (bskfit18.f) y
`axion_ns/axstar.py` resuelve TOV + Klein-Gordon (Ec. 2.6 del Paper I) para p_c = 100 MeV/fm^3 y
los modos (8a)-(8d) del Paper II con empalme a la Ec. (13) fuera (onda saliente o ligada):

- sin axión: 1.536 M_sun, 11.72 km; f, p1, p2 = 2.429, 7.171, 10.194 kHz (paper: 2.43, 7.15, 10.18);
- condensado entre 2.5 y 2.8 kHz; masas de la Tabla IV dentro del 0.35 %;
- Tablas I y II: desplazamientos nu/nu_0 - 1 dentro de 0.55 puntos (f: 0.33); mismos modos amortiguados
  salvo p1 a 8 kHz (en el umbral); |Im nu| dentro de un factor 2, tau_a mínimo 12.4 ms (paper: 12 ms).

La diferencia de fondo (≤0.3 %) viene de usar el ajuste analítico de BSk26 y no la tabla del grupo.
Cerca del centro las ecuaciones son rígidas: se integra con BDF hasta 0.5 km y luego con DOP853.

## Paso 10: con axión y métrica dinámica (etapa B, sin publicar)

`axion_ns/axgr.py` convierte fullgr.EXT a una matriz numérica (sin reescribir las ecuaciones a mano) y la
integra sobre el fondo de axstar: centro de la tercera parte, interior, cola del axión fuera (sistema
completo sin fluido) y, donde theta_0 ~ 1e-7, Zerilli + Ec. (13). El modo es un cero complejo de la
amplitud gravitacional entrante (Muller). Comprobaciones: sin axión, la conversión coincide con la del
paso 8 (1e-15) y el programa completo da el modo f del politropo de LMI (1.20162, Im al 0.2 %).

Resultados: sin axión f, p1, p2 = 1.961, 6.594, 9.872 kHz (Cowling: 2.429, 7.171, 10.194). Con axión, f
sube hasta un 7.7 % (Cowling 5.4 %) y su tau_GW baja de 0.174 a ~0.13 s; los p se mueven como en Cowling;
el máximo amortiguamiento (p2, 5 kHz, 12.8 Hz) coincide con Cowling, pero 1/tau_GW + 1/tau_a solo
acierta en orden de magnitud. Con delta theta apagado el sistema no es coherente: solo sus frecuencias.

Nota técnica: `simplify` de sympy no es determinista (depende de `PYTHONHASHSEED`) y, con algunas
semillas, su `trigsimp` no terminaba al separar el factor angular. Por eso `fullgr._divide` divide
como cociente de polinomios en (cos θ, sin θ) módulo cos² + sin² = 1. Probado con varias semillas.

## Estructura

| fichero | qué hace |
|---|---|
| `axion_ns/geometry.py` | Christoffel, Ricci, Einstein, divergencia covariante |
| `axion_ns/model.py` | acción, derivación de las ecuaciones de campo, fondo y linealización radial |
| `axion_ns/paper.py` | ecuaciones del Paper I transcritas |
| `axion_ns/nonradial.py` | linealización no radial en Cowling y ecuaciones del Paper II transcritas |
| `axion_ns/fullgr.py` | métrica perturbada, LD transcritas (cuatro fuentes) y el sistema extendido con axión |
| `axion_ns/fmodes.py` | integración numérica de las ecuaciones derivadas (modos de un politropo) |
| `axion_ns/bsk26.py` | ecuación de estado BSk26 (ajustes de Potekhin et al.) |
| `axion_ns/axstar.py` | estrella BSk26 con condensado axiónico y modos en Cowling con delta theta |
| `axion_ns/axgr.py` | sistema extendido convertido a números e integrado con axión (etapa B) |
| `axion_ns/compare.py` | comparación "módulo ecuaciones de fondo" |
| `axion_ns/step1…10_*.py` | los diez pasos (1–4: Paper I; 5–6: Paper II; 7–8: más allá de Cowling; 9–10: con axión) |
| `axion_ns/run_all.py` | ejecuta todo, resume y escribe `results/checks.json` |
| `report/` | genera la página HTML a partir de `results/checks.json` |
