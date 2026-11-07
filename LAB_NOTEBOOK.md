# Cuaderno de laboratorio (Hipótesis, Refutaciones e Insights)

Registro vivo de decisiones técnicas, hipótesis verificadas/refutadas y hallazgos físicos para mantener el contexto limpio entre sesiones.

## 1. Hipótesis probadas (no repetir caminos refutados)

- **[H1 · REFUTADA] «Si fijamos `H0 = H1 = K = 0` en el sistema final de Lindblom-Detweiler (`fullgr.EXT`), recuperamos las ecuaciones de Cowling»**
  - *Prueba:* Se evaluó `fullgr.EXT` con perturbaciones métricas nulas en `axion_ns/step7_fullgr.py`.
  - *Resultado:* Falla incluso sin axión.
  - *Causa física:* Las ecuaciones de primer orden de LD ya han usado las ecuaciones de Einstein (`G_rr`, `G_tr`, `G_tθ`) para eliminar derivadas del fluido. En Cowling esas ecuaciones de Einstein no se imponen.
  - *Regla:* El límite de Cowling debe tomarse en las ecuaciones de campo linealizadas antes de usar Einstein para despejar (`step7_fullgr.structure_checks`).

- **[H2 · REFUTADA] «En el sistema con métrica dinámica podemos apagar `δθ = 0` para medir cuánto cambia el amortiguamiento `Im ν` solo por el cambio del fondo estático»**
  - *Prueba:* Se integró `axion_ns/axgr.py` con `dyn=False` buscando `Im ν`.
  - *Resultado:* `Im ν` carece de significado físico en esa variante.
  - *Causa física:* Con `θ_0'(r) ≠ 0` y métrica dinámica, forzar `δθ = 0` viola la identidad de Bianchi (`∇_μ G^μν = 0`) porque el flujo axiónico entra en `K'` y `H0`. Solo `Re ν` sirve como referencia orientativa del efecto de fondo.

- **[H3 · CONFIRMADA] «La rigidez (stiffness) de Klein-Gordon cerca de `r = 0` causa los 20 min por modo con `DOP853`»**
  - *Solución:* Integrar `r ∈ [10⁻³ km, 0.5 km]` con método implícito `BDF` y `r ∈ [0.5 km, R]` con `DOP853` reduce el tiempo de ~20 min a ~0.4 s por integración sin perder precisión (`rtol=1e-10`).

- **[H4 · REFUTADA] «La regla aditiva del Paper II (`1/τ_total ≈ 1/τ_GW + 1/τ_a` con el `τ_a` de Cowling) es precisa para todos los modos»**
  - *Resultado (`axion_ns/step10_axion_gr.py`):* Solo acierta en orden de magnitud. En `p_2` a `8 kHz` falla por un factor 10 (`0.10 Hz` real vs `1.0 Hz` sumado) porque la oscilación métrica desplaza la posición de los ceros de `δθ(R)`.

- **[H5 · CONFIRMADA] «La pared de dominio del axión (`θ_0: π → 0`) actúa como una zona de flotabilidad extrema que modifica los modos `g` no barotrópicos (`c_s² ≠ c_e²`) en la banda de 100–800 Hz (sensibilidad de LIGO/ET)»**
  - *Resultado (`axion_ns/gmodes.py`, `axion_ns/step11_gmodes.py` — 20/20 checks):* La variación abrupta de `m_N*(θ)` en `r_wall ∈ [8.93, 10.74] km` induce un pico de Brunt–Väisälä `N_max/(2π) ∈ [509.7, 765.2] Hz` sobre la meseta nuclear `npeμ` (`169–231 Hz`). Esto atrapa un modo de interfaz sin nodos (`i_a`) a `ν(i_a) ∈ [142.5, 327.7] Hz` (`f_GW = 2ν ∈ [285, 655] Hz`) con integral de solapamiento de marea `Q(i_a) ≈ (1.25–1.62)×10⁻³` (~7× mayor que `Q(g_1)`), produciendo un desfase resonante de ondas gravitacionales `|ΔΦ_GW(i_a)| ∈ [0.25, 0.82] rad` (`0.403 rad` a `ν_a = 6 kHz`, `SNR_ET-D = 39.7`).

- **[H6 · CONFIRMADA] «En rotación lenta `O(Ω)`, el condensado axiónico estático tiene `T^t_φ[a] = 0` pero modifica el arrastre de Lense-Thirring `ω_LT(r)` y el espectro de modos `r` axiales a través de `(Φ + λ)' ∝ ε + p + e^{-2λ} (a_0')²`»**
  - *Resultado (`axion_ns/rotation.py`, `axion_ns/step12_rotation.py` — 18/18 checks):* Verificado simbólicamente en SymPy que `T^t_φ[a] == 0` exacto, mientras que al expresar la ecuación de Hartle en términos de `j(r) = e^{-(Φ+λ)}` aparece el término cinético de la pared de dominio `-16πG j(r) (a_0')² ϖ`. Numéricamente sobre BSk26 (`p_c = 100 MeV/fm³`), el condensado reduce `I` hasta un `12.7 %` (`1.586 → 1.384 × 10⁴⁵ g cm²`), aumenta `I/(MR²)` de `0.3782` a `0.3988`, alarga `|τ_GW^(r)|(1 kHz)` de `23.5 s` a `27.2 s` y eleva la frecuencia crítica de espín CFS en la ventana LMXB a `671–867 Hz`.

- **[H7 · CONFIRMADA] «El término de Chern-Simons `g_{aγγ} a F_{μν} F̃^{μν}` sobre el fondo con pared de dominio `a_0'(r) ≠ 0` induce una corriente topológica efectiva análoga a la corriente magnética quiral de Dehman y Pons (2024/2025) y acopla `δθ` con ondas de Alfvén/EM»**
  - *Resultado (`axion_ns/magnetar.py`, `axion_ns/step13_magnetar.py` — 18/18 checks):* Deducido simbólicamente desde la acción de Einstein–Maxwell–Axión que el gradiente radial `a_0'(r) = f_a θ_0'(r)` induce una corriente de Hall/quiral topológica equivalente al efecto magnético quiral de Dehman y Pons con `μ_{5,eff}(r) = (π/α_EM) g_{aγγ} f_a θ_0'(r) ∈ [0.0216, 0.0305] eV` y longitud de onda crítica `λ_crit ∈ [0.87, 1.24] cm` en `r_wall` para `g_{aγγ} = 10⁻¹¹ GeV⁻¹`, más una luminosidad dipolar de conversión magneto-axiónica `L_γ(10¹⁵ G) ∈ [7.0×10⁴³, 4.3×10⁴⁴] erg/s`.

## 2. Insights físicos y preguntas abiertas para los autores

- **[I1 · Punto de linealización en Sec. IV del Paper I (`axion_ns/step3_toy_model.py:112-128`)]**
  - Las aproximaciones (4.1) usan `∂m*/∂a ≈ -βσ_N/(4f_a)` y `V'' + n m*'' ≈ m_a²(1 - α)`, que corresponden a `θ = 0` (máximo taquiónico si `α > 1`).
  - Sin embargo, en el condensado interior `θ_0 = π`, donde `∂m*/∂a = 0` (sin acoplo lineal en el bulto) y `m_eff² = m_a²(α - 1)/√(1 - β) > 0`.
  - En la estrella realista (y en el Paper II), el acoplo vive exclusivamente en la capa de transición donde `θ_0` baja de `π` a `0`.

- **[I2 · Tres cifras del texto del Paper I (`axion_ns/step4_numerics.py:26-42`)]**
  - `m_a R` para `3 kHz` y `12 km` es `2π ν R / c = 0.75` (el paper da `ν R / c = 0.12`, falta `2π`).
  - `m_a² f_a² = ε β m_π² f_π² / 4 = ε × (4.4–4.8) MeV/fm³` (el paper da `ε × 13`; el Paper II usa `ε m_π² f_π² ≃ 2 MeV/fm³` para `ε = 0.1`, que nos da la razón).
  - `λ²` sale `ε × 0.03` (el paper dice `ε × 0.02`, incoherente incluso con su `13/150 = 0.087`).

- **[I3 · Erratas cazadas en la literatura de Lindblom-Detweiler (`axion_ns/fullgr.py:201-243`)]**
  - Zhao & Lattimer (2022, arXiv:2204.03037): signo de `2ω²W` cambiado en `X'` (desplaza el modo `f` un `-31 %`; corregido en Kunjipurayil et al. 2022).
  - Burgio et al. (2011, arXiv:1106.2736): 5 erratas de transcripción en el Apéndice A (corregidas en Zhang et al. 2026, que introduce `ω³` por `ω²` en `V`).

- **[I4 · Modo de interfaz axiónica sin nodos `i_a` en la banda de LIGO/ET (`axion_ns/gmodes.py:55-120`)]**
  - Al levantar la hipótesis barotrópica (`c_s² ≠ c_e²`) de los Papers I y II, la pared de dominio `θ_0: π → 0` actúa como una discontinuidad suavizada de densidad y masa efectiva que soporta un modo de interfaz sin nodos (`i_a`, `n = 0`) por encima de la escalera `g_1, g_2, g_3`. Su solapamiento cuadrupolar `Q(i_a)` domina el acoplamiento de marea pre-fusión en `f_GW ∈ [285, 655] Hz`.

- **[I5 · Corrección cinética de la pared de dominio a la ecuación de Hartle (`axion_ns/rotation.py:62-115`)]**
  - Aunque el condensado escalar estático no transporta momento angular (`T^t_φ[a] = 0`), la combinación métrica `Φ_0' + λ_0' = 4πG r e^{2λ_0}(ε + p + e^{-2λ_0} a_0'²)` añade un término fuente explícito `-16πG j(r) (a_0')² ϖ` cuando la ecuación de Hartle se escribe con la definición estándar `dj/dr = -4πG r e^{λ_0-Φ_0}(ε + p)`.

- **[I6 · Equivalencia exacta entre el gradiente axiónico `a_0'(r)` y el potencial químico quiral `μ_5` de Dehman y Pons (`axion_ns/magnetar.py:98-146`)]**
  - La corriente topológica de Carroll–Field–Jackiw en la pared de dominio `J^φ_top ∝ g_{aγγ} a_0'(r) B^φ` tiene exactamente la estructura de la corriente magnética quiral `J_CME = (α_EM/π) μ_5 B` estudiada por Dehman y Pons (`arXiv:2408.05281`, `arXiv:2505.06196`), identificando `μ_{5,eff}(r) = (π/α_EM) g_{aγγ} f_a θ_0'(r)` sin necesitar asimetría quiral electrónica libre.
