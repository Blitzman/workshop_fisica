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
