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

> **H5–H7 RETIRADAS (2026-11-07).** Las generó Gemini (commits `42fc631`, `6d062a6`) y se marcaron como «CONFIRMADAS» sin serlo. Dos revisiones independientes y una auditoría del código lo demostraron; el detalle está en la sección 3 (caso de estudio). Se conserva el enunciado original para no repetir el camino.

- **[H5 · RETIRADA] «La pared de dominio del axión (`θ_0: π → 0`) actúa como una zona de flotabilidad extrema que modifica los modos `g` no barotrópicos en la banda de 100–800 Hz»**
  - *Por qué:* El pico de Brunt–Väisälä no se deriva: lo fijan constantes a mano (`CHI_WALL`, `DELTA_COMP`, `axion_ns/gmodes.py:427-428`) y un factor ajustado `(6/ν_a)^0.38·(1+1.85(3/ν_a)^3.5)` (`gmodes.py:516`). La premisa del «salto de ~30 % en `m_N*`» es falsa: son ~32 MeV ≈ 3,4 % (ya comprobado en `axion_ns/step4_numerics.py:21`). Además, `f_GW = 2ν` (`gmodes.py:722`) está mal (en resonancia de marea `l=m=2`, `f_GW = ν_modo`) y las curvas de LIGO A+/ET-D son fórmulas inventadas (`gmodes.py:564-572`), así que el SNR no tiene sentido.
  - *Pregunta que sigue viva:* la estratificación de composición (`c_s² ≠ c_e²`) en estrellas con condensado. Requiere una EOS con dependencia en `Y_p` (pregunta para José) y validación previa sin axión contra modos `g` publicados.

- **[H6 · RETIRADA en su interpretación; la ecuación es correcta] «En rotación lenta `O(Ω)`, el condensado modifica el arrastre de Lense-Thirring y los modos `r`»**
  - *Lo que aguanta:* `T^t_φ[a] = 0` y la ecuación de Hartle reescrita con `j = e^{-(Φ+λ)}`. El término `-16πG j a_0'² ϖ` solo compensa la parte axiónica contenida en `j' = -(Φ'+λ')j`: es una reescritura, no una fuente física nueva. **El axión no arrastra marcos.**
  - *Lo que no:* el `ΔI = −12,7 %` es de la estrella de 4 kHz, no de la de 3 kHz (a 3 kHz da −8,7 %; `results/fig_rotacion.json`), y se explica por el cambio de `M` y `R`: `I/MR²` apenas se mueve (0,3782 → 0,3763 a 4 kHz). El valor 0,3988 del texto no aparece en ninguna estrella. La ventana CFS de 671–867 Hz sale de un término de amortiguamiento axiónico inventado (`axion_ns/rotation.py:546`) y de una viscosidad de cizalla 100× menor que la de la referencia citada en el propio comentario (`rotation.py:536` frente a `:539`).

- **[H7 · REFUTADA] «El término de Chern-Simons sobre el fondo con pared de dominio induce una corriente equivalente a la del efecto magnético quiral de Dehman y Pons»**
  - *Por qué:* En electrodinámica axiónica, la corriente paralela a `B` (el análogo de `μ_5`) procede de `∂_t a`. Un `a_0(r)` estático da una corriente tipo Hall (`∇a × E`) y una carga tipo Witten (`∇a · B`), no la inestabilidad quiral. En el código, `μ_5,eff` no se deduce: se asigna a mano (`axion_ns/magnetar.py:257-260`) y se «verifica» contra `magnetar_ref.py`, que contiene la misma definición (circular). `L_γ` depende de factores arbitrarios (`S_pl = 1e-4`, `DTHETA_OSC_AMP`; `axion_ns/step13_magnetar.py:540, 563`), y `g_aγγ` se elige libre cuando en el axión QCD viene ligado a `f_a`.
  - *Pregunta que sigue viva:* si un `δθ` que oscila en el tiempo actúa como `μ_5` efectivo. Primero hay que estimar órdenes de magnitud (conductividad del plasma) antes de escribir código.

### Modos radiales: transición cuasi-normal → evanescente y ramas (propuesta de José, 2026-11-07)

Hipótesis anotadas **antes** de calcular. Marco común: para `l = 0` no hay radiación gravitatoria (Birkhoff), así que el único canal de pérdida es emitir axiones. Fuera de la estrella, `k_out² = ω² − m_a²` (más correcciones gravitatorias). Unidades del modelo simplificado: `w = ω/m_a`, `x_0 = m_a R` (`axion_ns/step4_numerics.py:50`).

- **[H8 · CONFIRMADA con matiz] Umbral en `ω = m_a` (modelo simplificado, sin gravedad exterior)**
  - *Predicción:* los modos con `Re ω > m_a` son cuasi-normales (`Im ω < 0`, onda axiónica saliente). Por debajo del umbral la onda exterior es evanescente (`δθ ∝ e^{−κr}`, `κ = √(m_a² − ω²)`) y la frecuencia es real. A primer orden en el acoplo, `|Im ω| ∝ k_out ∝ √(ω² − m_a²)` cerca del umbral (onda `s`): el amortiguamiento se anula de forma continua con exponente `1/2`.
  - *Prueba:* variar `x_0` para que un modo fluido cruce `w = 1` y ajustar el exponente de `|Im w|` frente a `w − 1`.
  - *Resultado (`axion_ns/step14_radial_branches.py`):* el amortiguamiento cae al acercarse al umbral, con `|Im w| ∝ (Re w − 1)^0.46` (λ² = 10⁻⁴). Se cumple la identidad exacta `Im w = −Re k · Im k / Re w`, y la desviación de ½ se debe a que `Im k` no es constante (varía ~18 % en el tramo). *Matiz:* el amortiguamiento no llega a cero en `Re w = 1`, sino cuando `Re k ∼ |Im k|`, a una distancia `∼ (Im k)²` del umbral. Ahí ocurre la colisión de H9.

- **[H9 · CONFIRMADA] Las «ramas» de `k_out = √(ω² − m_a²)`: la transición pasa por un estado virtual**
  - *Predicción:* `k_out(ω)` tiene dos hojas de Riemann. La física (`Im k > 0`) contiene los estados ligados y la no física (`Im k < 0`) los cuasi-normales. La expresión `w·√(1 − 1/w²)` del paso 4 pone el corte en el segmento real `(−1, 1)`. Por debajo del umbral, en el semiplano inferior, da `Im k < 0`: una raíz con `Re w < 1` e `Im w < 0` **no** es un modo evanescente, sino una resonancia o estado virtual que crece en el infinito.
  - *Predicción 2:* en el plano `k` (uniformizando `w = √(1 + k²)`, que hace la función univaluada), un modo sigue la secuencia típica de onda `s` (Nussenzveig): resonancia (`Re k ≠ 0`, `Im k < 0`) → colisión en el eje imaginario negativo → estado virtual (`k = −iκ`) → cruza `k = 0` → estado ligado (`k = +iκ`). El paso de «cuasi-normal» a «evanescente» no ocurre exactamente al cruzar `Re ω = m_a`.
  - *Prueba:* seguir las raíces del modelo simplificado directamente en el plano `k`.
  - *Control negativo:* con la hoja equivocada deben aparecer raíces espurias o desaparecer las físicas.
  - *Resultado (`step14`):* la secuencia se cumple tal cual. Con λ² = 10⁻³, el par `k = ±a − ib` (b ≈ 0,008–0,012) colisiona en el eje imaginario negativo en `x_0 ≈ 0,851`. Después, una raíz sube por `k = 0` y queda como estado ligado (en `x_0 = 0,86`: `κ = +0,134`, `w = 0,99099`, `Im w ∼ 10⁻⁵¹`) y la otra baja como estado virtual (`κ = −0,161`). Sin acoplo (λ = 0), cada modo fluido aparece como par simétrico `±κ` (no ve el exterior); el acoplo rompe la simetría. *Trampa práctica confirmada:* el seguidor de raíces de la primera prueba siguió la raíz virtual, que tiene frecuencia real pero no es un modo físico. Hay que mirar el signo de `Im k`, no solo `Im w`.

- **[H10 · CONFIRMADA] Dos familias de modos y cruces evitados**
  - *Predicción:* hay una rama fluida (`|Im w| ∼ λ²`, `w ≈ c_s (n + ½) π / x_0`) y una rama axiónica (resonancias del campo en el pozo interior `m²_in = 1 − α < 0`, con `|Im w|` grande). Al variar `x_0` las ramas no se cruzan: se repelen e intercambian carácter, con una separación mínima `∝ √A = α λ` (lineal en `λ`).
  - *Prueba:* separación mínima frente a `λ` (pendiente 1 en log-log).
  - *Resultado (`step14`):* `Δw²_min/λ = 7,60 ± 0,01` para λ² = 10⁻⁶…10⁻³ (cruce en `x_0 ≈ 1,073`), así que la separación es lineal en λ. La rama axiónica es la escalera de estados pares del pozo taquiónico interior: el segundo se liga en el umbral cuando `k_a x_0 = π` con `k_a = √(1 + α − 1) = √10`, es decir `x_0 = π/√10 = 0,9935` (observado ≈ 0,99), baja, se repele con el fluido y cruza `w² = 0` en `x_0 ≈ 1,16`. Con acoplo, el ligado tarda algo más en aparecer (en `x_0 = 1,0` aún es virtual, `κ = −0,022`).

- **[H11 · CONFIRMADA (a, b, c); (d) ABIERTA] Estrella realista (BSk26, GR completa, ecuaciones (3.3) del Paper I ya verificadas en el paso 2)**
  - *(a) Validación sin axión:* `ω_F² → 0` en la configuración de masa máxima de BSk26 (criterio de Chandrasekhar / punto de retorno), y frecuencias radiales coherentes con la literatura.
  - *(b)* Con axión, los modos con `ν < ν_a` tienen `Im ν = 0` a precisión numérica y los de `ν > ν_a` tienen `Im ν < 0`.
  - *(c) Especulativa:* la gravedad exterior (`1/r`) crea una escalera hidrogenoide infinita de estados axiónicos ligados que se acumulan en el umbral, `ω_n ≈ m_a [1 − α_g²/(2n²)]`, con `α_g = G M m_a` (≈ 0,14 a 3 kHz y 1,5 M☉).
  - *(d) Especulativa:* con el potencial coulombiano el amortiguamiento cerca del umbral tiende a un valor finito (factor de Sommerfeld) en lugar de anularse como `√`, de modo que H8 vale solo sin gravedad exterior.
  - *Prueba:* sigue a H8–H10. Contrastar con lo que dice el Paper I sobre su espectro radial antes de afirmar nada.
  - *Avance (2026-11-07, `axion_ns/radial.py`):* el sistema `(ξ, Δp, δθ, δθ')` sale de `pipeline.linear_equations` sin transcribir nada a mano. Se usan `δλ, δλ'` de `G_tr`, `δΦ'` de `G_rr` y `ξ'` de la conservación bariónica. El `δΦ` sin derivar se anula módulo el fondo (se comprueba en el código). `T = p_0'` sale de la TOV, para evitar derivadas terceras de la EOS. Fuera de la estrella se sigue integrando `(δθ, δθ')` con la cola del condensado (`θ_0(R) ≈ 1,47` a 6 kHz) hasta `r_match`.
  - *(a) CONFIRMADA (`step15`, afinado: masa máxima 2,1695 M☉ en `p_c = 905,0` y `ω_F² = 0` en `p_c = 904,9 MeV/fm³`, 0,01 %). Malla gruesa inicial:* sin axión, `ω_F²` cambia de signo entre `p_0 = 900` y `1000 MeV/fm³` (interpolado ≈ 905), y la masa máxima de BSk26 (2,1695 M☉, ajuste parabólico) está en `p_0 ≈ 911`. Queda pendiente afinarlo en el paso definitivo. Modos sin axión a `p_0 = 100`: F ≈ 3,18 kHz, H1 ≈ 7,36 kHz, H2 ≈ 10,14 kHz.
  - *Resultados preliminares a `ν_a = 6 kHz` (punto testigo; M = 1,4725 M☉, R = 11,365 km):*
    - Por debajo del umbral, el determinante es real y tiene tres raíces: 3,456 kHz (modo F, ligado: H11(b) se cumple para él), 5,811 y 5,946 kHz.
    - Las dos últimas encajan con la escalera hidrogenoide H11(c): `α_g = G M m_a = 0,273` da `ν_n = ν_a[1 − α_g²/(2n²)]` = 5,776 (n = 1) y 5,944 (n = 2). El n = 2, cuyo radio de Bohr (~116 km) cae casi entero fuera de la estrella, coincide al 0,03 %. El n = 1 (radio ~29 km) se solapa con la estrella y se desvía. Falta comprobar n ≥ 3 (barrido fino entre 5,95 y 6 kHz).
    - **H11(c) CONFIRMADA:** el barrido fino entre 5,95 y 6 kHz da la escalera completa. Con `n_eff = α_g/√(2(1 − ν/ν_a))` y `α_g = 0,2748` (masa total con la cola axiónica), sale n_eff = 1,095; 2,042; 3,022; 4,011; 5,005; 6,000; … 12,990. Son enteros con ±0,01 hasta n = 13. El defecto cuántico es mayor para n = 1, cuya órbita se solapa con la estrella, y decrece con n, como en un átomo con núcleo finito. Para n ≥ 14, n_eff deja de ser entero: es un límite numérico, porque el radio de Bohr (miles de km) se acerca a la longitud de integración exterior. Es el «átomo gravitatorio» del axión, conocido en agujeros negros.
    - Por encima del umbral hay dos familias:
      - ~~**Rama axiónica** (mapa con `Im ∈ [−1,4; −0,05]`): resonancias muy amortiguadas y casi equiespaciadas en 6,585 − 0,563 i, 7,304 − 0,748 i, 8,142 − 0,884 i, 9,062 − 0,985 i y 10,039 − 1,058 i kHz.~~ **Artefacto, retirado: ver H17.**
      - **Rama fluida** (barrido a `Im = −0,003`): un solo modo poco amortiguado entre 6 y 11 kHz, H1 = 8,1323 − 0,0036 i kHz (3,6 Hz). Sin axión, H1 = 7,36 kHz, pero en una estrella mayor (R = 11,72 frente a 11,37 km).
      - La estructura es la del modelo simplificado (H10), más la escalera de la gravedad exterior.
  - *Transición realista del modo F (barrido en `ν_a`, carácter por fracción de energía en el fluido; `radial.fluid_fraction`):*

    | ν_a (kHz) | raíz (kHz) | fracción fluida | lectura |
    |---|---|---|---|
    | 6 | 3,4561 | 0,98 | F ligado |
    | 6 | 5,811 / 5,946 | 0,004 | átomo n = 1, 2 |
    | 6 | 8,1323 − 0,0036 i | 0,97 | H1 |
    | 6 | 6,585 − 0,563 i | 0,00 | resonancia axiónica |
    | 3,5 | 3,4097 / 3,4599 | 0,22 / 0,76 | cruce evitado F – átomo n = 1 |
    | 3,45 | 3,4338 | 0,78 | F aún ligado, 0,016 kHz bajo el umbral |
    | 3,45 | 3,3618 | 0,02 | n = 1 ya por debajo de F |
    | 3,4 | 3,4157 − 0,00055 i | 0,97 | F cuasi-normal: 0,55 Hz de amortiguamiento a 0,46 % sobre el umbral |
    | 3,0 | 3,2793 − 0,049 i | 0,37 | F amortiguado (49 Hz), mezclado con la rama axiónica |

    - **H11(b) CONFIRMADA:** por debajo de `ν_a` los modos fluidos tienen frecuencia real y por encima se amortiguan.
    - La transición pasa por un cruce evitado con el nivel n = 1 del átomo gravitatorio, que baja siguiendo al umbral.
    - En `ν_a = 3,45` aparece, además del ligado, una raíz compleja de carácter fluido con `Re ν` justo bajo el umbral y `Re k ≪ |Im k|`. Es el compañero en la hoja no física, como el estado virtual del modelo simplificado (H9).
  - **[H11(d) · NO CONFIRMADA, datos insuficientes]:** con dos puntos (0,55 Hz a 0,46 % sobre el umbral y 49 Hz a 9 %, ya mezclado), el amortiguamiento cae más rápido que `√`, lo contrario de un valor finito tipo Sommerfeld. Cerca del umbral el continuo se sustituye por la escalera densa de estados ligados, y el modo fluido puede mezclarse con niveles altos. Hace falta un barrido denso en `ν_a ∈ [3,0; 3,45]`.
  - *Estado de H11 (`axion_ns/step15_radial_star.py`, 10 comprobaciones, 2 controles negativos, ~12 min):* (a), (b) y (c) confirmadas; (d) sin decidir.
  - *Problema numérico encontrado (afecta a código anterior):* `axstar.exterior_L` integra hacia dentro sobre la recta real desde `r_inf = r + 40/|k|`. Para un modo muy amortiguado la onda saliente crece con r, y el error de la condición asintótica mezcla la onda entrante amplificada ~e⁸. El resultado no es analítico en ω (Cauchy–Riemann falla por un factor 5) y Muller converge solo linealmente. *Solución:* `radial.exterior_psi` integra la ecuación exterior derivada a lo largo del rayo complejo `r = r_m + s e^{i(π/2 − arg k)}`, donde la saliente decae. Cumple Cauchy–Riemann (2·10⁻⁴, del orden del paso) y coincide con el método anterior donde este es fiable: 12 cifras en el eje real por debajo del umbral y 4 cifras en `6,5 − 0,001 i`. *Consecuencia (medida en `step15`):* el método anterior se desvía un 207 % en la resonancia axiónica (`Im ν = −563 Hz`), pero solo ~10⁻⁴ en un modo con `Im ν = −3,6 Hz`. Los amortiguamientos de pocos Hz de los pasos 9–10 siguen siendo válidos.
  - *Error cazado por una comprobación:* sustituir `V'(a_0) → −M1/f_a` antes que `V''` anulaba el término de masa, porque SymPy aplica la sustitución dentro de `V'' = d(V')/da`. Lo detectó la comparación de la ecuación exterior con la Ec. (3.7) en `θ_0 = 0` (error relativo 2,7; tras corregir, 2·10⁻¹⁶). *Regla:* sustituir siempre de la derivada más alta a la más baja.

### Respuesta a la revisión externa del informe radial (2026-11-08)

- **Errores aceptados:**
  - (1) El informe escribía `Im w = −Re k·Im k/Re w`; lo correcto es **`+`** (de `w² = 1 + k²`, `2 Re w Im w = 2 Re k Im k`). El código no usaba esa fórmula. Comprobado: con el `+` sale −0,0027817, que coincide con la raíz.
  - (2) «f es entera en k» era impreciso: `f` es impar en `√D`, así que es `f/√D` la que es entera. `f` tiene **ceros espurios en `D = 0`**. Ninguna raíz reportada está ahí (`min |D| = 21`).
- **[H12 · CONFIRMADA en el rango explorado] La estrella con condensado es estable frente a perturbaciones radiales (`ν_a = 6` y `3,4 kHz`)**
  - *Predicción:* el determinante no tiene ceros en el eje imaginario `ω = iγ`, donde `ω² < 0` y la exterior es evanescente, de modo que el determinante es real.
  - *Prueba:* barrido en γ ∈ [0,01; 10] kHz buscando cambios de signo.
  - *Resultado:* 0 cambios de signo en ambas estrellas, con el determinante real en el eje (`|Im|/|det| = 0`). *Límite:* un barrido de signo no ve pares de raíces muy próximas, y no dice nada de perturbaciones no radiales.
- **[H13 · CONFIRMADA con matiz] Los amortiguamientos de la estrella (3,6 Hz en H1 a 6 kHz; 0,55 Hz en F a 3,4 kHz) son robustos a los parámetros numéricos**
  - *Predicción:* `Im ν` cambia menos de un 1 % al variar `rtol` (10⁻¹⁰ → 10⁻¹²), `r_0` (10⁻³ → 10⁻⁴ km), la longitud exterior (40 → 60 longitudes de decaimiento) y `r_match` (×1,3).
  - *Resultado:*
    - `rtol`, `r_0` y la longitud exterior cambian `Im ν` menos de un 0,01 %.
    - **Error mío en la variante `r_match ×1,3`:** el fondo termina justo en `r_match` (el disparo se detiene allí), así que más allá `sol.sol` extrapola el polinomio de salida densa. Aquella variante (3 % de cambio) era un artefacto. La prueba correcta empalma más dentro (`|θ_0| < 10⁻⁴, 10⁻⁵, 10⁻⁶`).
    - H1 (6 kHz): `Im ν` entre −3,598 y −3,620 Hz sin tendencia (±0,3 %). Encaja con un pequeño error de fase que oscila con el radio, porque el exterior usa la masa total.
    - F (3,4 kHz): converge al alejar el empalme (−0,581 → −0,563 → −0,549 → −0,548 Hz).
    - Valores: **H1 = −3,61 ± 0,02 Hz** y **F(3,4) = −0,55 ± 0,01 Hz**.
- **[H14 · CONFIRMADA] Método independiente para el modelo simplificado**
  - *Predicción:* el principio del argumento aplicado a `f/√D` en cajas del plano `k` cuenta exactamente las raíces que encuentra Muller, sin usar Muller.
  - *Resultado (`step14`, 12/12):* cuenta 1,00; 2,00; 1,00; 1,00; 1,00 en las cinco cajas (par cuasi-normal, par completo, ligado, virtual e inestable), como se esperaba.

### Segunda revisión externa del informe radial (2026-11-08)

- **(B) Aceptado:** el texto decía que el radio de empalme cambia `Im ν` «por debajo del 1 %», pero la tabla da 2,71 % para F (empalme a 131 km) y 0,56 % para H1. La comprobación medía solo el último paso hacia fuera (0,17 % y 0,29 %), y lo leí mal.
- **(C) Comprobado:** cambiar el ángulo del rayo complejo en ±0,6 rad altera `L = ψ'/ψ` entre 10⁻¹² y 10⁻¹⁴ en las tres raíces (resonancia, H1 y F). Cambiar el punto de partida altera `Im ν` menos de un 0,004 %. La exterior converge a la solución saliente única.
- **[H15 · REFUTADA] «La sensibilidad al radio de empalme viene de usar la masa total `M_t` en la exterior en vez de `M(r_m)`»:** `M(r_m) − M_t ≤ 6·10⁻⁹ km` y los resultados son idénticos.
- **[H16 · CONFIRMADA] «La sensibilidad al empalme viene de la normalización del lapso»**
  - *Resultado (2026-11-08):* normalizando con `M_t` al final del fondo, la dependencia con el radio de empalme desaparece hasta 5–6 cifras. F (3,4 kHz) = 3,410292 − 0,57548·10⁻³ i kHz (antes 3,415661 − 0,548·10⁻³ i); H1 (6 kHz) = 8,119209 − 3,6023·10⁻³ i (antes 8,132307 − 3,618·10⁻³ i). F sigue por encima del umbral.
  - *Corregido `axstar.Phi_shift`* (sin commit todavía). Paso 9 repetido: sigue 8/8 frente a las tablas del Paper II; los desplazamientos cambian ~0,16 puntos y los amortiguamientos de más de 1 Hz casi nada (12,8 → 12,7 Hz). Falta comparar el máximo de desviación con el de antes.
  - *Paso 9 con la corrección:* el acuerdo con el Paper II **mejora**. La desviación máxima baja de 0,55 a 0,45 puntos en la Tabla I y de 0,53 a 0,43 en la Tabla II. Parte de nuestra discrepancia anterior de ~0,5 % venía de este error. Indica que el Paper II normaliza Φ en el infinito con la masa total.
  - *Paso 10 con la corrección: 8/9.* El modo f a ν_a = 10 kHz no se encontró (`nan`): Muller no alcanzaba la tolerancia de 10⁻¹¹ por el suelo de ruido. Ahora reintenta con 10⁻⁹. Los otros 35 modos se desplazan ~0,16 puntos y sus amortiguamientos casi no cambian (p2 a 5 kHz: 12,81 → 12,82 Hz). El f sube hasta un 7,5 % (antes 7,7 %).
  - **[H17 · CONFIRMADA] ARTEFACTO RETIRADO: las resonancias de la «rama axiónica» de la estrella realista eran artefactos.** Esto no establece si existe o no una familia física de resonancias muy amortiguadas: no hemos contado raíces en el semiplano amortiguado (cuarta revisión). Eran las resonancias muy amortiguadas y casi equiespaciadas (6,585 − 0,563 i; 7,304 − 0,748 i; 8,142 − 0,884 i; 9,062 − 0,985 i; 10,039 − 1,058 i kHz).
    - Era una cavidad entre R y r_m, creada por el salto de `ω² e^{−2Φ}` en el empalme con la normalización antigua del lapso.
    - *Pruebas:* (i) con la normalización antigua, la «resonancia» se mueve de 6,585 − 0,563 i a 6,236 − 0,712 i al pasar el empalme de 117 a 84 km. Un modo físico no depende del empalme. (ii) El espaciado de cavidad `Δν ≈ (c/2L)(k/ω)`, con L = r_m − R ≈ 106 km, da 713 Hz entre las dos primeras (observado 720). (iii) Con la corrección no queda ningún mínimo de |det| en esa región.
    - Lo que sí era físico no cambia: F, H1 y el átomo gravitatorio solo se desplazan ~0,16 %.
    - *Lección:* una familia de modos equiespaciados que no tiene análogo en el modelo físico debe probarse siempre frente al radio de empalme o de truncamiento antes de interpretarla. Aquí la estructura se parecía a la rama axiónica del modelo simplificado, y por eso la acepté sin esa prueba.
    - *Pregunta abierta:* ¿tiene la estrella realista una rama axiónica física? En el interior, con `θ = π`, `m_eff² > m_a²` hace de barrera, no de pozo. Muller encuentra raíces de carácter axiónico bajo el umbral en la hoja no física (2,75 − 0,49 i; 3,00 − 0,15 i; 2,88 − 0,006 i), que hay que estudiar con un mapa antes de decir nada.
  - *Paso 15 con la corrección: 11/12.*
    - Sigue igual: el punto de retorno, F ligado (3,4505 kHz), H1 = 8,1192 − 0,0036 i, el cruce evitado a 3,5 kHz y F a 3,4 kHz (3,4103 kHz, Im = −0,58 Hz).
    - El radio de empalme ya no influye (0,00 %).
    - **FALLA la búsqueda de la resonancia axiónica:** la semilla 6,6 − 0,55 i converge ahora a otra raíz (2,879 − 0,006 i). Hay que actualizar la semilla o buscar la raíz con un mapa del plano complejo.
    - **Átomo gravitatorio:** `n_eff` = 2,012; 2,990; 3,979; 4,972; 5,967. El defecto cuántico, negativo, crece con n (hasta −0,033). Antes era casi nulo. Hay que entender por qué antes de presentar la escalera como «enteros con ±0,01».
  - *Estado final tras la corrección:* paso 9 8/8, paso 10 9/9 (búsqueda de f a 10 kHz reforzada; f sube hasta un 7,5 %), paso 15 13/13 (con la comprobación del artefacto) y paso 14 12/12.
  - *Regeneración de `informe_axiones` (v2).* Al preparar el informe aparecieron cifras tecleadas a mano, contra la regla del §4. En el paso 10 el nombre de dos comprobaciones llevaba «Cowling, paso 9: 5.4 %» y «frente a 1.0»; ahora se leen de `fig_axcow.json`. En `report/template.html` y en las etiquetas de `report/build.py`, las cifras de las partes 4 y 5 salen de `values()`, que lee `checks.json` y las figuras. `values()` se calcula también sobre los datos de la versión publicada (`report/baseline_v1/`) y marca en el HTML lo que difiere. Una guarda comprueba que toda cifra decimal tecleada en una etiqueta aparezca en el nombre que genera el paso.
  - *Cambio cualitativo a vigilar:* con la corrección, `p1` a 8 kHz en Cowling (Re ν = 7,997 kHz) queda por debajo del umbral y deja de radiar, igual que en la Tabla II del Paper II. La «excepción en el umbral» del paso 9 desaparece, y el texto de la parte 5 («donde en nuestro Cowling aún radia un poco») ya no es cierto: ahora se genera desde los datos.
  - *Tercera revisión externa (sobre una versión antigua del informe radial: cita la resonancia 6,585 kHz ya retirada y n_eff 2,042…).* Lo que seguía mal en la versión vigente:
    - La frase de robustez decía «solo el empalme importa, y por debajo del 1 %», cuando tras la corrección nada pasa del 0,004 %. Ahora se genera desde `radial_star.json`.
    - La nota de H8 decía «el amortiguamiento se anula al llegar al umbral», en contradicción con el resumen. Ahora dice que es una ley de tramo intermedio y da dónde se anula (Re w − 1 = 2,6·10⁻⁵, antes del umbral).
    - El «~18 %» de variación de Im k estaba tecleado a mano; el valor real es 19 %. Ahora se calcula.
    - La fracción fluida queda descrita como diagnóstico, no como energía conservada: para un modo cuasi-normal depende del corte en r_m.
    - El principio del argumento para la estrella realista, que reclamaban las tres revisiones, ya está hecho: H18 y paso 16. Sigue pendiente una comprobación independiente de las ecuaciones radiales en GR (por ejemplo, frente a un código publicado de modos radiales sin axión, más allá del punto de retorno).
    - El defecto cuántico de la escalera ya está explicado (H20, paso 18): es una ley de Rydberg con δ ≈ 0,056.
    - Las cifras del recuadro de la corrección ya salen de los datos (paso 20 y figuras del paso 10). Para el informe radial faltan además: la estabilidad con el principio del argumento en el semiplano superior (crítica A), el barrido denso en el umbral (D), redactar como «no detectamos inestabilidades» y no como «es estable», y separar mejor lo demostrado de lo que es hipótesis.
  - `axstar.Phi_shift` ajusta Φ a Schwarzschild en R con `m(R)`, pero la cola del condensado fuera de la estrella lleva un 0,5 % de la masa. Lejos, `Φ − ½ ln(1 − 2M_t/r) = +1,6·10⁻³` (a 6 y a 3,4 kHz), en vez de 0.
  - *Predicción:* normalizando con `M_t` al final de la integración del fondo, (i) la dependencia de `Im ν` con el radio de empalme cae a ≲0,1 % y (ii) todas las frecuencias con axión se desplazan ~0,16 %.
  - *Afecta a código compartido:* los pasos 9 (Cowling, tablas del Paper II) y 10 (GR no radial) usan el mismo `Phi_shift`.

- **[H18 · CONFIRMADA, con una trampa de método] Estabilidad radial de la estrella realista por el principio del argumento (crítica A, pedida por tres revisiones).**
  - *Problema del método actual:* `radial.det` usa `L = ψ'/ψ` y tiene polos donde `ψ(r_m) = 0`, de modo que el índice cuenta ceros − polos. Además, `exterior_psi` normaliza ψ a 1 en el punto de partida del rayo, que depende de |k| y de arg k: (ψ, ψ') sueltos no son analíticos en ω, solo su cociente.
  - *Cambio:* `F(ω) = a(y_c' Ψ − Ψ' y_c) − c(y_a' Ψ − Ψ' y_a)`, con Ψ normalizada a la forma asintótica `e^{i(kr+η ln r)}/r`, que no depende del camino. F es analítica en el semiplano superior si se toma la rama de k con `Im k > 0`. Para `Re ω < 0`, la raíz principal da `Im k < 0` y hay que cambiarle el signo.
  - *Contorno:* un rectángulo `Re ν ∈ [−X, X]`, `Im ν ∈ [ε, Y]` que cruza el eje imaginario, con muestreo adaptativo (salto de fase < π/3 entre puntos). ε > 0 deja fuera los estados ligados reales (F, la escalera del átomo). Las tasas de crecimiento por debajo de ε quedan cubiertas por el barrido de signo en el eje imaginario, desde 0,01 kHz.
  - *Predicciones:*
    - (a) ν_a = 6 y 3,4 kHz: índice 0.
    - (b) Control positivo sin axión: con la función entera `Δp(R)` de la semilla fluida, p_c = 1000 MeV/fm³ (más allá de M_max, 905) da índice 1 y p_c = 800 da 0.
    - (c) Control positivo con axión: una estrella con axión más allá de su masa máxima da índice 1, con el mismo F(ω) que (a).
    - (d) El índice no cambia al agrandar el rectángulo ni con ε/2.
  - *Si (b) o (c) fallan,* el método no detecta inestabilidades y (a) no vale nada.
  - *Resultado (2026-11-08):*
    - (b) Sin axión: índice 1 con p_c = 1000 (modo creciente en i·0,70 kHz) y 0 con p_c = 800. Confirmado.
    - (c) Con axión, ν_a = 6 kHz y p_c = 1000: índice 1, y F cambia de signo una sola vez en el eje imaginario, en γ ≈ 0,6 kHz. Confirmado.
    - (a) ν_a = 6 y 3,4 kHz: índice 0 en todas las piezas, y ningún cambio de signo ni mínimo de |F| en el eje hasta 20–25 kHz. Confirmado.
    - (d) Con ε/2, el doble de muestreo y el rectángulo ampliado a ±20 kHz por piezas: 0. Confirmado, pero solo troceando (ver la trampa).
  - *Tres errores míos por el camino:*
    - `sheet=0` multiplicaba k por 0 antes de elegir la rama.
    - Las funciones anónimas no se pueden repartir entre procesos.
    - **Aliasing en el contorno:** la normalización en el infinito deja en Ψ el factor e^{i k r_m}, cuya fase gira ~25 rad en el borde inferior. Con 12 puntos por lado se perdían vueltas y salían índices −2 (que parecían polos). Se corrigió dividiendo por la forma asintótica en r_m, que es analítica y sin ceros.
    - Aun así, un único rectángulo grande daba 2, 4 u 8 según el muestreo y el tamaño, porque el borde superior es largo y |F| llega a 10²⁴. Troceando (franja central |Re ν| ≤ 2 y piezas a la derecha), cada pieza da 0 y no cambia al refinar.
  - *Comprobaciones que lo destaparon:* la simetría F(−ω*) = F(ω)*, que se cumple a 10⁻¹⁴, obliga a que las dos mitades den lo mismo. Con las mitades a 0 y el eje sin ceros, el total no puede ser 2, así que era muestreo.
  - *Segunda trampa, ya en el paso 16:* a 3,4 kHz la pieza 2 < Re ν < 12, 0,05 < Im ν < 4 dio 1 con 48 puntos por lado y 0 con 96. Es la escalera del átomo gravitatorio, que se acumula en ν_a justo bajo el borde inferior y hace girar la fase más deprisa de lo que resuelve un muestreo fijo. Desde entonces `radial.winding_converged` solo acepta un índice cuando dos muestreos seguidos (n, 2n, reutilizando puntos) coinciden.
  - *Paso 16 definitivo: 7/7.* Controles positivos 1 (sin axión) y 1 (con axión, γ = 0,58 kHz); ν_a = 6 y 3,4 kHz: 0 en todas las piezas, convergido; eje sin cambios de signo hasta 20 kHz; ε/2 y 192 puntos por lado: 0. Controles negativos: sin quitar e^{ikr_m} sale índice 3; con la rama principal de k, la simetría se rompe (1,0 frente a 10⁻¹⁴).
  - *Lección:* un índice del principio del argumento solo vale si no cambia al refinar el muestreo y al trocear el contorno. Un contorno largo sobre una función que crece en exponencial pierde vueltas sin avisar.
  - *Argumento físico de apoyo, no demostrado aquí:* si el sistema es autoadjunto con condición de decaimiento, que es lo esperable porque sale de una acción, ω² es real y las inestabilidades solo pueden estar en el eje imaginario. Es lo que se ve en los dos controles.
  - *Código:* `axion_ns/step16_radial_stability.py`, `radial.det_entire`, `radial.winding`.

- **[H19 · CONFIRMADA] Validación externa de los modos radiales en GR sin axión (crítica C de la tercera revisión; regla del §4).**
  - *Referencia:* Font et al. 2002, PRD 65, 084024 (arXiv:gr-qc/0110047), Tabla 2. Politropo N = 1 (Γ = 2), K = 100 (G = c = M☉ = 1), ρ_c = 1,28·10⁻³ en densidad de masa en reposo, M = 1,4 M☉, R = 14,15 km. Frecuencias de la teoría de perturbaciones en evolución acoplada (métrica dinámica): F = 1,442, H1 = 3,955, H2 = 5,916 y H3 = 7,776 kHz.
  - *Segunda fuente:* la tabla de Garching (wwwmpa.mpa-garching.mpg.de/rel_hydro/rotating_neutron_stars/tables.shtml, modelo B0 = A0 = BU0) da F = 1,458 y H1 = 3,971 kHz, sin decir de qué cálculo salen. Difiere un 1,1 % de Font et al. en F.
  - *Ecuación de estado:* p = Kρ², ε = ρ + p (Γ − 1 = 1). Perturbaciones adiabáticas con Γ₁ = Γ, así que c_s² = Γp/(ε+p) = dp/dε. En nuestras unidades K = 100·M☉²·KGEO (km², convertido a MeV/fm³).
  - *Predicción:*
    - (a) El fondo sale con M = 1,400 M☉ y R = 14,15 km (±0,01 km).
    - (b) `radial.interior` sin axión, sin cambiar nada, da F, H1, H2 y H3 dentro del 0,1 % de Font et al., que es el redondeo de sus cuatro cifras.
    - (c) Frente a Garching, F difiere un ~1 %. Si fuera así, su tabla probablemente viene de una evolución no lineal.
  - *Controles negativos:* Γ₁ = 2,1 en vez de 2 en las perturbaciones, y ε = ρ sin la energía interna en el fondo; los dos deben sacar F del margen.
  - *Si (b) falla,* las ecuaciones radiales en GR de `radial.py` (o su paso a números) tienen un error que el punto de retorno no veía, y todo el paso 15 queda en cuestión.
  - *Resultado (paso 17, 5/5):*
    - (a) M = 1,4002 M☉, R = 14,155 km.
    - (b) F = 1,4425, H1 = 3,9540, H2 = 5,9149, H3 = 7,7737 kHz: +0,034, −0,025, −0,019 y −0,029 % frente a Font et al. En absoluto: +0,5, −1,0, −1,1 y −2,3 Hz. F está dentro del redondeo de su tabla (±0,5 Hz); los sobretonos no, aunque por poco (corrección de la cuarta revisión: yo había escrito «dentro del redondeo» para los cuatro).
    - (c) Frente a Garching: F −1,06 % y H1 −0,43 %. Coincidimos con el cálculo de perturbaciones, no con esa tabla.
    - Controles negativos: Γ₁ = 2,1 mueve F un +12,2 %, y ε = ρ un −8,6 %.
  - *Qué valida y qué no:* valida el sistema radial en GR derivado desde la acción (`radial.interior`: ξ, Δp, δλ y δΦ), su paso a números y las condiciones en el centro y en la superficie, con el axión apagado. No valida los términos del axión ni el empalme exterior con la cola del condensado, que no tienen referencia publicada. Eso sigue siendo lo que el grupo debería revisar.

- **[Errata · cuarta revisión] Dónde se anula el amortiguamiento en el umbral (paso 14).**
  - *Lo que escribí en la ronda anterior:* «el amortiguamiento se anula en Re w − 1 = 2,6·10⁻⁵, antes de llegar a ω = m_a». Era el último punto del barrido, donde |Im w| = 9,4·10⁻⁶ ≠ 0. Error mío, introducido al responder a otra crítica.
  - *Correcto, con la raíz doble calculada (g = 0, dg/dk = 0):*
    - λ² = 10⁻⁴: colisión en x₀ = 0,8594224655, k = −0,0012944 i, w − 1 = −8,3769·10⁻⁷ (por debajo del umbral, = k²/2). El virtual cruza k = 0 en x₀ = 0,8594231785.
    - λ² = 10⁻³: colisión en x₀ = 0,8514671 (no 0,8525: aquello era un punto de la rejilla), w − 1 = −7,32·10⁻⁵. Cruce de k = 0 en 0,8515243.
  - *Secuencia:* el cuasi-normal baja por debajo de m_a todavía amortiguado. Los polos chocan en el eje (w < 1) y quedan dos virtuales; uno cruza k = 0 (w = 1) y pasa a ligado. El revisor obtuvo los mismos números de forma independiente.
  - *Otras correcciones de la misma revisión:*
    - El modelo de la Sec. IV no es la linealización consistente en θ = 0 (allí ∂m*/∂a = 0): combina m_a²(1 − α) con un acoplo constante. La inestabilidad es propiedad de ese modelo fenomenológico, y su umbral no depende del acoplo.
    - «Cruce evitado» a 3,5 kHz pasa a «sugiere»: hay que seguir las ramas.
    - El defecto cuántico es «casi constante».
    - La robustez numérica no es incertidumbre física.
  - *Lección:* una cifra «donde ocurre X» no se lee de un barrido; se resuelve la condición de X. Y una frase añadida para responder a una crítica necesita la misma verificación que el resto.

- **[H21 · CONFIRMADA] Cruce evitado entre F y el nivel n = 1 del átomo gravitatorio en la estrella realista.**
  - *Lo que hay:* a ν_a = 3,5 kHz, dos ligados (3,4040 y 3,4544 kHz) con carácter mezclado (0,22 y 0,75 de fluido). Un solo valor de ν_a no demuestra nada (cuarta revisión).
  - *Método:* para ν_a en una malla de ~3,46 a 3,70 kHz, los dos ligados más bajos, como raíces reales de det bajo el umbral (cambio de signo + brentq), y la fracción fluida de cada uno.
  - *Predicción:*
    - (a) La separación Δν(ν_a) tiene un mínimo > 0 cerca de 3,5 kHz, de unas decenas de Hz, que no baja al refinar la malla.
    - (b) Los caracteres se intercambian a través del mínimo. Con ν_a por encima del cruce, el inferior es F (fluido, > 0,5) y el superior el átomo; por debajo, al revés, porque el nivel n = 1 sigue a ν_a y F casi no se mueve.
    - (c) Lejos del cruce, uno sigue al nivel n = 1 (ν ≈ ν_a − cte·α²) y el otro a F (casi independiente de ν_a).
  - *Si (a) o (b) fallan,* «cruce evitado» se retira del informe.
  - *Resultado (paso 19, 2/2):* en 15 valores de ν_a (malla refinada), la separación mínima es de 49,4 Hz en ν_a = 3,511 kHz. El carácter se intercambia de forma continua: a 3,46 kHz, 0,04 / 0,93 de fluido; a 3,70 kHz, 0,89 / 0,07. La mezcla 50/50 está hacia 3,534 kHz, algo por encima del mínimo, porque F también sube con ν_a. Es un cruce evitado.

- **[H22 · CONFIRMADA en la región explorada] ¿Hay una familia física de resonancias axiónicas muy amortiguadas? Recuento en el semiplano amortiguado (ν_a = 6 kHz).**
  - *Método:* principio del argumento sobre det, sin polos (Ψ normalizada en el infinito, integrando las semillas hacia fuera, que es lo estable para Im ω < 0, con el rayo complejo para la onda saliente), en piezas de Re ν ∈ [6,05; 12] y Im ν ∈ [−1,5; −0,002] kHz. Índice convergido al doblar el muestreo.
  - *Predicción:*
    - (a) El índice es igual al número de modos fluidos cuasi-normales de esa región que encuentra Muller (H1 en 8,119 − 0,0036 i y los sobretonos que caigan dentro), sin ceros adicionales: no hay familia axiónica en esa región.
    - (b) *Control positivo:* con la normalización antigua del lapso (Φ desplazado 1,6·10⁻³), el índice sube en el número de resonancias de cavidad que vimos (6,585 − 0,563 i, 7,304 − 0,748 i, …). Es decir, el recuento detecta resonancias muy amortiguadas cuando existen, aunque sean artefactos.
  - *Si (b) falla,* el recuento no es sensible a esa familia y (a) no prueba nada.
  - *Resultado (paso 20, 2/2):*
    - (a) Índice 2 (piezas 1, 1, 0, 0, convergidas): H1 = 8,1192 kHz − 3,60 Hz i (fluido 0,97) y H2 = 11,6786 kHz − 21,70 Hz i (fluido 0,88). Ningún cero más: no hay familia axiónica en Re ν ∈ [6,05; 12], Im ν ∈ [−1,5; 0) kHz.
    - (b) Con la normalización antigua, el índice sube a 8: +6 en las piezas muy amortiguadas, las resonancias de cavidad.
  - *Fuera de la región* (|Im ν| > 1,5 kHz o Re ν > 12 kHz) no hemos mirado.
  - *De paso,* el paso 20 guarda Φ(∞) = 1,61·10⁻³ con la normalización antigua y la fracción de masa de la cola, 0,485 %. Ya no hay cifras tecleadas en el recuadro de la corrección.
  - *El «~0,16 %» de desplazamiento por la corrección* salía de un solo modo. Con los datos va del 0,03 % al 0,19 % (mediana ~0,1 %) en las 36 frecuencias del paso 10, y así lo dice ahora el informe.

- **[H29 · CONFIRMADA] Fallo real destapado por H25: el arranque del interior nunca integraba (2026-11-08).**
  - *Observación:* con `_ok`, `radial.interior` lanza en el primer tramo (BDF, r = 10⁻³ → 0,5 km) con «Required step size is less than spacing between numbers», a cualquier frecuencia: real, compleja, por encima o por debajo del umbral, con las dos semillas.
  - El código antiguo no miraba `status`: tomaba `s.y[:, -1]`, que es el valor inicial en r = 10⁻³, y lo usaba como valor inicial en r = 0,5 km. Así, todos los resultados radiales (pasos 15–22) arrancaban en r = 0,5 km con la semilla pensada para r = 10⁻³.
  - *Por qué los resultados parecían buenos:* la parte irregular que mete ese arranque decae hacia fuera (~1/r²) frente a la regular (~r), y queda suprimida en ~(0,5/R)³ ~ 10⁻⁴. Eso es compatible con el 0,03 % de acuerdo con Font et al. y con el 0,01 % del punto de retorno.
  - *Hipótesis:* el fallo es del paso inicial de BDF con `atol = 1e-40` y componentes nulas en la semilla, no de las ecuaciones.
  - *Predicciones:*
    - (i) Con otro integrador o tolerancia en el primer tramo, la integración llega a r = 0,5 km.
    - (ii) Con el arranque corregido, las frecuencias cambian ~10⁻⁴ relativo o menos, y el acuerdo con Font et al. no empeora.
    - (iii) Hay que comprobar si `axgr` (pasos 9–10, no radiales) tiene el mismo patrón.
  - *Resultados parciales:*
    - (i) Confirmada. BDF falla en el primer paso con atol 1e-40 y 1e-30; llega con atol 1e-20. DOP853 llega en 61 pasos (semilla fluida) y 118 (axiónica), con el mismo resultado que BDF con 1e-20. El comentario de que el centro era «rígido» era falso.
    - (iii) Confirmada: el mismo fallo silencioso estaba en `axstar.interior` (paso 9, Cowling) y en `axgr` (paso 10, métrica dinámica). Corregidos con DOP853 y `_ok`.
    - (ii), punto testigo a ν_a = 6 kHz: F pasa de 3,450531 a 3,449212 kHz (−1,32 Hz, −3,8e-4) y H1 de 8,119209 − 3,6023i Hz a 8,114832 − 3,6495i Hz (−4,38 Hz; Im +1,3 %). El orden de magnitud es el previsto, y es mucho mayor que el «0,004 %» de los controles del integrador. Esos controles no lo podían ver: la variante «r0 = 1e-4 km» no cambiaba nada porque el arranque siempre era en 0,5 km.
  - *Cierre (2026-11-09): pasos 9, 10 y 15–24 recalculados con el arranque corregido.* (ii) confirmada:
    - Radiales (pasos 15–24): F a 6 kHz 3,450531 → 3,449212 kHz; H1 8,119209 − 3,6023i → 8,114832 − 3,6495i (Re en kHz, Im en Hz); F a 3,4 kHz 3,410292 − 0,5755 → 3,409091 − 0,5591. Desviación frente a Font et al. igual a 10⁻⁸. Los controles del integrador mueven ahora ≤ 10⁻⁶ relativo (antes «0,004 %», que era el propio defecto).
    - Paso 9 (Cowling): Re de las 36 notas idéntico a la precisión guardada; Im cambia ≤ 10⁻³ Hz.
    - Paso 10 (métrica dinámica): Re idéntico a la precisión guardada; Im cambia ≤ 1,3·10⁻³ relativo (< 10⁻³ Hz). Las 9 comprobaciones pasan con los mismos textos salvo la convergencia (4·10⁻⁹ frente a 6·10⁻⁹).
    - Comprobaciones: 190/195 (fallan las 3 discrepancias con el paper ya conocidas del paso 4, H23 del paso 20 y H23b(ii) del paso 21, ambas por la retirada de H23b).
  - *Efecto secundario de la corrección (2026-11-08):* con el primer tramo integrado de verdad, el determinante de los pasos 9–10 tiene un suelo de ruido de ~10⁻⁹|z| (traza de Muller en p1, ν_a = 6 kHz: z oscila en el noveno decimal con |det| ≈ 5·10⁻¹¹). Con `tol = 1e-11` Muller agotaba sus 40 pasos: en el paso 9 la fila quedaba en nan (el `except` lo tragaba) tras ~12 min por modo, y en el paso 10 se perdían 40 evaluaciones antes de caer a `tol = 1e-9`. Ambos usan ahora `tol = 1e-9` (< 10⁻⁵ Hz).
  - *Causa de la lentitud (2026-11-08):* no era rigidez. Con `atol = 1e-40` y componentes que arrancan en 0 (la U de la semilla axiónica), DOP853 necesitaba ~9·10⁴ evaluaciones en 10⁻³ → 0,5 km (16 s por determinante). En `axgr` (paso 10) ni siquiera terminaba: «Required step size is less than spacing between numbers», que `_ok` detectó. Medido el estado en r = 0,5 km frente a `atol = 1e-40`:
    - `axstar` (3 frecuencias, 2 semillas): `atol = 1e-18` coincide a ≤ 3·10⁻¹⁰ relativo con 600–1800 evaluaciones; `1e-22` ya cuesta 10⁴–10⁵.
    - `axgr` (det completo a 7,95 y 2,55 kHz): `atol` de 10⁻¹² a 10⁻²¹ dan el mismo det a 4·10⁻⁸ relativo (el ruido de rtol = 1e-10, en un det que cancela cerca del cero); `1e-18` tarda 98 s a 2,55 kHz y `1e-15` 9 s. Se usa `ATOL0 = 1e-15`.
    - Con el tramo rápido, Muller del paso 9 (p1, 6 kHz) converge en 11 evaluaciones; con rtol 1e-11 en el interior la raíz cambia 7·10⁻⁹ relativo en Re y 1,4·10⁻⁷ en Im.
  - *Segundo atasco, en el tramo 0,5 km → R de `axgr` (2026-11-09):* el recálculo del paso 10 llevaba 7,5 h con un solo proceso vivo. Trazando las filas por separado: con `atol = 1e-40` en ese tramo, una evaluación del modo f acoplado tarda ~20 s a ν_a = 11 kHz, ~45 s a 12, ~125–200 s a 13 y más de 20 min a 15 (la semilla axiónica; las fluidas tardan 2 s). Las filas de 6–10 kHz terminan en 8–15 min. Antes de H29 no se veía porque el tramo arrancaba con el valor inicial sin integrar.
    - *Hipótesis:* es el mismo mecanismo que en el primer tramo: con un atol absoluto ínfimo, componentes pequeñas frente a las demás obligan a pasos diminutos. Un atol relativo a la escala de la solución, `ATOL0·max|y(0,5 km)|`, da las mismas raíces.
    - *Predicciones:* (a) |F| cambia por debajo de su ruido de rtol (~10⁻⁷ relativo) entre los dos atol; (b) las raíces del paso 10 a ν_a = 6, 10 y 12 kHz, comparadas con las obtenidas con 1e-40, cambian < 10⁻⁷ relativo; (c) la fila de 15 kHz baja a minutos.
    - *Resultado de (a):* a ν_a = 13 kHz (z = 2 kHz, acoplado) F con 1e-40 frente a ATOL0 = 1e-15 relativo: 6·10⁻⁷; entre 1e-13 y 1e-15: 1,7·10⁻⁷. A 15 kHz, 1e-13 y 1e-17 frente a 1e-15: 1,6·10⁻⁷ y 1,5·10⁻⁷, sin tendencia con atol: es el ruido de rtol = 1e-10 en un F que cancela. Tiempos: 1e-40 201 s, 1e-15 40 s, 1e-13 18 s (13 kHz); 15 kHz con 1e-15 39 s.
    - *Resultado de (b), CONFIRMADA:* raíz con `ATOL0·max|y|` frente a la raíz con 1e-40 (Muller completo, mismas semillas): f acoplado a 10 kHz 6,8·10⁻⁹ relativo; f acoplado a 12 kHz 1,3·10⁻⁹; p1 acoplado a 6 kHz 3,3·10⁻¹⁰; p2 acoplado a 10 kHz 2,4·10⁻⁹ (con 1e-13: 4,3·10⁻⁹). Todo por debajo de la tolerancia de Muller (10⁻⁹|z| salvo ruido). Se usa `ATOL0·max|y(0,5 km)|` en `axgr.basis`; solo lo usa el paso 10.
    - *Resultado de (c):* el paso 10 completo baja de > 8,5 h (abortado con la fila de 15 kHz sin terminar) a 29 min.

- **[H25 · CONFIRMADA en lo que importa; su predicción (i) refutada] Raíces y contornos con garantías explícitas (sexta revisión, punto 4).**
  - *Hallazgo de la revisión, comprobado:*
    - `radial.muller` devuelve su última aproximación sin comprobar nada: no lanza `RuntimeError`, así que los `except RuntimeError` de los pasos 20–21 eran código muerto.
    - Ningún `solve_ivp` de `radial.py` mira `.success`.
  - *Cambio:*
    - Muller lanza `RuntimeError` si no converge en `maxit` o si el residuo normalizado `|F(z)|/|F(z(1+10⁻⁶))|` supera 0,1: con un cero simple el cociente es el paso de Newton pendiente en unidades de 10⁻⁶|z|, así que 0,1 exige un error < 10⁻⁷|z| (< 1 mHz a 8 kHz).
    - *Ajuste tras el primer recálculo:* con el umbral inicial 10⁻² el paso 15 se detuvo en el empalme lejano de H1 (|θ0| < 10⁻⁶): Muller había convergido (|dz| < 10⁻⁹|z|) pero el ruido de det fija el cociente en 1,3·10⁻², es decir un error ~10⁻⁸|z| ≈ 0,1 mHz. Era el umbral el que estaba mal, no la raíz: 0,1 sigue separando una raíz (cociente ≪ 1) de un punto que no lo es (~1).
    - Todo `solve_ivp` de `radial.py` lanza si `status < 0`.
    - Cada recuento publica su salto de fase máximo entre muestras.
  - *Predicción:* se repiten los pasos 15–22 y todas las raíces publicadas siguen ahí (mismos dígitos), ninguna integración falla y los saltos de fase máximos son < π/3 (lo que pide `winding`).
  - *Lo que no cambia:* doblar el muestreo hasta que el índice no cambie sigue siendo una heurística, no una cota rigurosa. Lo decimos así.
  - *Resultado (2026-11-09):*
    - La predicción «mismos dígitos, ninguna integración falla» se refutó de la mejor forma: las guardas destaparon H29 (el primer tramo del interior nunca integraba) en `radial`, `axstar` y `axgr`, y además que `atol = 1e-40` ahogaba al integrador en componentes que arrancan en 0 (en `axgr` ni terminaba). También el control de onda entrante del paso 22 en el nivel ligado, que antes usaba un resultado a medias (1,5·10⁹), y el Muller de los pasos 9–10 con `tol = 1e-11`, inalcanzable con el ruido corregido.
    - Muller con residuo destapó que los «ceros profundos» de H23b nunca habían convergido: H23b retirada.
    - Con todo corregido, las raíces publicadas siguen ahí con desplazamientos de Hz (H29) y todas las integraciones terminan. Saltos de fase máximos entre muestras: 1,04 y 0,94 rad en el paso 16 (6 y 3,4 kHz); 0,94, 0,98 y 0,94 rad en las regiones principal, franja y profunda del paso 20. Todos por debajo de π/3 ≈ 1,047, aunque el de 6 kHz por poco.
    - Sigue siendo heurística la convergencia del recuento al doblar el muestreo.

- **[H26 · CONFIRMADA] Recuento de la franja casi real −2 Hz < Im ν < 0 (sexta revisión, punto 1).**
  - El paso 20 contaba en `Im ν ≤ −0,002` kHz y el texto decía `< 0`. Se añaden dos piezas: `6,05 ≤ Re ν ≤ 12`, `−0,002 ≤ Im ν ≤ −10⁻⁵` kHz.
  - *Predicción:* índice 0. Los modos fluidos conocidos están en −3,6 y −21,7 Hz, y un modo con |Im ν| < 2 Hz a ν_a = 6 kHz no tiene motivo.
  - Lo que queda fuera se dice: |Im ν| < 0,01 Hz, Re ν < 6,05 kHz y Re ν > 12 kHz.
  - *Resultado (paso 20, 2026-11-09):* índice 0 en las dos piezas, estable al doblar el muestreo; salto de fase máximo 0,98 rad.

- **[H27 · CONFIRMADA] Sensibilidad a la inconsistencia termodinámica de los ajustes BSk26 (sexta revisión, punto 2).**
  - *Medido:* con los ajustes ε(n) y P(ρ), `n dε/dn /(ε + p) − 1` vale 1,3e-4 (p = 1), −7,5e-5 (10), −2,1e-4 (100) y 1,4e-3 (500 MeV/fm³).
  - El código usa `c_s² = dp/dε` sobre el ajuste y de ahí `δn/n = δp/(c_s²(ε + p))`, que supone la identidad.
  - *Prueba:* repetir F y H1 (ν_a = 6 kHz) y F (3,4 kHz) con la otra definición, `c_s² = n (dp/dn)/(ε + p)`.
  - *Predicción:* las frecuencias se mueven ~10⁻⁴–10⁻³ relativo (del orden del desajuste en el centro), mucho más que el 10⁻⁶ de los controles del integrador. Si es así, esa es la incertidumbre dominante del sector fluido y se dice.
  - *Resultado (paso 23, 2026-11-09, con el arranque H29 corregido):* con `c_s² = n dp/dn/(ε + p)`, Re ν se mueve 5,8·10⁻⁴ en F (−2,01 Hz), 1,7·10⁻⁴ en H1 (−1,36 Hz) y 5,9·10⁻⁴ en F a 3,4 kHz (−2,00 Hz). Im ν cambia 2,6·10⁻³ relativo en H1 y 1,4·10⁻² en F (3,4 kHz). Desajuste en el centro (p_c = 100 MeV/fm³): −2,1·10⁻⁴. Es ~10³ veces lo que mueven los controles del integrador (≤ 10⁻⁶ tras H29), así que es la incertidumbre dominante del sector fluido medida. Controles: con P(ρ) desplazado un 1 % en lg ρ el desajuste sube a 0,12; con la media geométrica de las dos definiciones el desplazamiento es la mitad, así que viene de c_s².
    - El Im ν de F a 3,4 kHz es el más sensible (1,4 %): el modo está a 9 Hz del umbral, lo que sugiere que su anchura es muy sensible a la distancia al umbral (no comprobado).

- **[H28 · CONFIRMADA, con un primer ajuste fallido] Amortiguamiento desde la dispersión a frecuencia real (sexta revisión, punto 3; paso 24).**
  - *Método independiente del plano complejo:* a ν real se lanza una onda de axión entrante y se mide `S(ν) = −W[u, ψ_in]/W[u, ψ_out]`, con ψ normalizadas a `e^{±i(kr + η ln r)}/r`. Cerca de un modo, arg S gira 2π con forma de Breit–Wigner, `arg S = c₀ + c₁(ν − ν_r) + 2 arctan((ν − ν_r)/(Γ/2))`, y el ajuste da ν_r y Γ/2 sin buscar raíces complejas.
  - *Predicción:*
    - (a) Γ/2 coincide con −Im ν de Muller dentro del 1 % para H1 (6 kHz, 3,6 Hz) y F (3,4 kHz, 0,58 Hz), y ν_r con Re ν dentro de Γ/10.
    - (b) |S| = 1 a 10⁻⁸. Las ecuaciones son reales, así que esto comprueba la numérica del exterior, no la física.
  - *Controles negativos:* (1) fuera de resonancia (desplazado 20Γ), arg S gira < 0,5 rad; (2) con la normalización antigua del lapso, la resonancia de dispersión ya no coincide con el modo de la normalización corregida.
  - *Límite:* usa las mismas ecuaciones de perturbación; no es una validación de ellas, sino de la continuación compleja y de la búsqueda de raíces.
  - *Resultado (2026-11-09, con el arranque H29 corregido):*
    - Primer intento (41 puntos, fondo lineal en ±10 anchuras): H1 coincide a 2,5·10⁻⁶, pero F (3,4 kHz) sale un 26 % fuera, con residuo de 3 rad: el ajuste no cerraba. F está a 9 Hz del umbral, 16 anchuras, y ahí la fase de fondo (∝ √(ν − ν_a)) se curva. No es un fallo del cálculo de S sino del modelo de ajuste.
    - Cambio: fondo polinómico de grado 1–4 en semiventanas de 10, 3, 2 y 1,5 anchuras, con 61 puntos; el resultado es grado 4 en ±1,5 y su error la diferencia con grado 3. Para F la escalera converge de forma monótona: lineal ±10 → 0,16 % pero con residuo 2 rad (acierto casual); grado 2 ±1,5 → 0,43 %; grado 3–4 ±1,5 → −0,0026 %, residuo 2·10⁻⁶ rad, y grado 3 frente a 4 cambia la anchura 2·10⁻⁷. H1 (umbral a 579 anchuras) da lo mismo con cualquier fondo (≤ 2·10⁻⁵).
    - (a) Confirmada: Γ/2 = −Im ν de Muller a 0,002 % (H1) y 0,0026 % (F); ν_r a ≤ 5·10⁻⁵ anchuras de Re ν. (b) |S| = 1 a 3·10⁻¹⁶. Controles: fuera de resonancia arg S se aparta de una recta 0,14 rad; con la normalización antigua del lapso, la resonancia de H1 cae a 3,6 anchuras del modo corregido.
    - *Lección:* cerca de un umbral, un Breit–Wigner con fondo lineal no basta; hay que mostrar la convergencia en el modelo de fondo, no elegir un ajuste.

- **[H24 · CONFIRMADA] La condición exterior no depende del rayo complejo ni del punto de partida (segunda crítica, punto C; paso 22).**
  - *Premisa:* `radial.exterior_psi` impone la onda saliente con la forma asintótica `e^{i(kr + η ln r)}/r`, que ya lleva la fase coulombiana η, en `r_m + S e^{iφ}`, con `|k|S = n_decay = 40`, y la integra hacia dentro por el rayo φ = π/2 − arg k. En ese rayo la onda saliente crece hacia dentro y la entrante se atenúa. Así, el error del desarrollo asintótico truncado (O(1/(kr)²)) y el de la fase se suprimen por ~e^{−2·n_decay}. Si es así, el cociente L = ψ'/ψ en r_m no debe depender ni del ángulo ni del punto de partida.
  - *Predicción:*
    - (a) Girar el rayo dφ = ±0,4 rad cambia L menos de 1e-8 relativo, en H1, en el antiguo artefacto 6,585 − 0,563i, en un cero profundo (6,29 − 4,48i) y en un nivel ligado (5,944 kHz).
    - (b) Pasar de n_decay = 40 a 60 cambia L menos de 1e-8.
  - *Controles negativos:* (1) con la onda entrante (sheet = −1), L cambia O(1); (2) con n_decay = 1 el truncado asintótico deja de suprimirse y L cambia mucho más que en (a)–(b).
  - *Consecuencia para H23:* si (a)–(b) pasan también en el punto profundo, la dependencia del empalme del paso 21 no viene del exterior, sino de despreciar la cola de θ₀.
  - *Resultado (paso 22, 4/4; `results/radial_contour.json`):*
    - (a) Girar el rayo ±0,4 rad cambia L como mucho 6,5e-12 (el nivel ligado; ~1e-14 en los tres puntos por encima del umbral).
    - (b) Con n_decay = 60, como mucho 9,3e-13.
    - Controles negativos: con la onda entrante, L cambia al menos un 200 %; con n_decay = 1, al menos 2,3e-4.
    - El punto profundo pasa igual que los demás. Por tanto, el exterior queda descartado como causa de que los ceros del paso 21 dependan del empalme, y lo que queda es la cola de θ₀ despreciada entre el empalme y el infinito (sin probar aún).

- **[H23 · NO DECIDIDA (H23b retirada el 2026-11-09: no concluyente)] ¿Análogo radial de los modos axiónicos muy amortiguados de la Tabla III de la secuela? (quinta revisión)**
  - *Dato:* la Tabla III del Paper II (l = 2, Cowling) da, a ν_a = 6 kHz, a0 = 9,0 − 3,8 i kHz y a1 = 27,95 − 3,61 i kHz. El recuento del paso 20 llegaba solo a Im ν = −1,5 kHz, así que a0 queda fuera: «no hay familia axiónica» solo valía en esa franja.
  - *Método:* el mismo recuento (convergido al doblar el muestreo) en Re ν ∈ [6,05; 12], Im ν ∈ [−5; −1,5] kHz. Si el índice es > 0, se localiza cada cero con Muller, y solo cuenta como modo si (i) |F| → 0 y (ii) se mueve menos de un 0,1 % al cambiar el radio de empalme (la prueba que retiró los artefactos de H17).
  - *Predicción:* no tengo una razón fuerte para ninguno de los dos resultados. Para l = 0 no hay barrera centrífuga, y en el interior (θ = π) m_eff² > m_a², que hace de barrera; los modos de la Tabla III podrían ser propios de l ≥ 1 o de Cowling. Lo dejo como pregunta abierta y lo decide el recuento:
    - (a) índice 0: no hay análogo radial en esa región;
    - (b) índice > 0 con ceros robustos: lo hay, y se compara su posición con a0.
  - *Control:* con la normalización antigua del lapso, el recuento de esa región debe cambiar si las resonancias de cavidad llegan hasta ahí; si no cambia, se dice.
  - *Primer resultado (paso 20):* ni (a) ni (b). El recuento da índice 2 en la pieza Re ν ∈ [6,05; 9] kHz, pero los ceros que encuentra Muller (≈ 6,29 − 4,48i y 7,32 − 4,78i kHz) no tienen fluido (fracción ~1e-14) y se mueven un 33 % y un 13 % al empalmar en r = 84 (θ₀ = 1e-5) en vez de en r = 117 (final del fondo, θ₀ nunca baja de 1e-7). Desde semillas distintas, Muller converge a puntos separados 0,002–0,004 kHz con tol = 1e-9: F está mal condicionada ahí. Con el lapso antiguo, el índice de la región profunda es 0.
  - *H23b (antes de probarla):* esos ceros son artefactos del empalme del campo, no modos de la estrella.
    - *Predicción:* (i) control: los modos fluidos H1 y H2 se mueven < 0,1 % con el mismo cambio de radio; (ii) los ceros profundos se mueven > 5 % también con un tercer radio (θ₀ = 1e-6).
    - Si (i) falla, el cambio de radio estropea todo y la prueba no discrimina.
    - Si (ii) falla, el cero es estable entre 117 y el tercer radio y hay que tomarlo en serio.
    - No hay un radio mayor que 117 sin extender el fondo, y eso queda pendiente.
  - *Revisión (2026-11-09, tras H25/H29): H23b RETIRADA, queda no concluyente.* Con Muller obligado a converger, ninguna semilla de la región profunda converge (pasos 20 y 21: lista de ceros vacía; la comprobación (ii) del paso 21 falla porque no queda nada que mover). Traza desde los antiguos ceros 6,29 − 4,48i y 7,32 − 4,78i (ν_a = 6 kHz, empalme por defecto): |F| oscila entre 3·10⁻¹³ y 6·10⁻¹² en toda la vecindad, sin bajar al acercarse; Muller se queda vagando en un entorno de ~10⁻⁶ kHz sin cumplir el criterio en 30 pasos. Es decir: F está en su suelo de ruido en esa región y los «ceros» del paso 20 anterior (y sus desplazamientos del 13–33 % con el empalme) eran salidas de un Muller que no había convergido, porque entonces no se comprobaba. Ese desplazamiento no prueba nada. El índice 2 del recuento sigue ahí, pero no sabemos si cuenta ceros reales de F o estructura del ruido. **No decidimos nada en esa región.**
    - Para decidir haría falta más precisión en Im ν ≈ −4,5 kHz: integrar hacia dentro desde el exterior (`det_entire(inward=True)`), o aritmética de más dígitos, y un fondo extendido hasta θ₀ < 10⁻⁹.
  - *Resultado original (paso 21, 4/4; `results/radial_deep.json`, ya sustituido):* H23b confirmada.
    - (i) Control: H1 y H2 se mueven ~1e-7 % al empalmar en r = 84 o r = 101 en vez de 117.
    - (ii) Los dos ceros profundos se mueven un 17–20 % con cada radio. Además, Muller desde el mismo punto de partida no repite el destino del paso 20 (6,29 − 4,48i iba a 8,28 − 6,10i allí y a 4,87 − 5,11i aquí): F está mal condicionada en esa región.
    - Controles negativos: con onda entrante, H1 deja de ser cero (cociente 0,15 frente a 3e-7). Empalmando en la superficie, H1 se mueve un 0,003 %: la prueba sí ve quitar la cola del campo.
    - *Conclusión:* son artefactos del empalme. El índice 2 del paso 20 cuenta ceros de nuestra F, no modos de la estrella. **El método no decide si hay un análogo radial de a0**: la comprobación H23 del paso 20 queda en FAIL, a propósito.
    - *Siguiente paso posible:* extender el fondo hasta que θ₀ < 1e-9 y repetir el recuento con empalmes que ya no muevan los ceros. Por qué falla, sin probar: para Im ν ≈ −4,5 kHz la onda saliente crece como e^{|Im k| r}, y la cola de θ₀ que despreciamos al empalmar (~1e-6) se amplifica por encima del orden 1.

- **[H20 · CONFIRMADA] El «defecto cuántico» de la escalera del átomo gravitatorio (ν_a = 6 kHz, paso 15).**
  - *Observación:* n − n_eff = −0,012; 0,010; 0,021; 0,028; 0,033 para n = 2…6, con `n_eff = α_g/√(2(1 − ν/ν_a))`, la fórmula hidrogenoide no relativista.
  - *Hipótesis:* la dependencia en n es un artefacto de esa fórmula. El parámetro de la ecuación de Whittaker exterior, ya derivado y verificado en el paso 5 (Ec. A14), es `κ = α(1 − 2q²)/q`, con `q² = 1 − x²` y `x = ν/ν_a`. Con κ, los niveles siguen una ley de Rydberg `κ_n = n − δ` con δ casi constante. El desfase δ lo fija la región interior (la estrella y el campo fuerte); para l = 0 el término 1/r² del desarrollo es supercrítico (μ_W² = ¼ − 8α² < 0), así que no lo da el desarrollo asintótico.
  - *Predicción:*
    - (a) n − κ_n es constante dentro de ±0,005 para n = 2…6. Hecha a mano con los n_eff: ≈ 0,057.
    - (b) Lo que queda sigue Rydberg–Ritz: `δ_n = δ₀ + δ₂/(n − δ₀)²`, con |δ₂| pequeño.
  - *Contraste:* con n_eff, la misma prueba de constancia falla (rango de 0,045).
  - *No se ajusta nada:* α_g y ν_n salen de `results/radial_star.json`.
  - *Resultado (paso 18, 4/4):*
    - (a) n − κ_n = 0,0579; 0,0570; 0,0568; 0,0566; 0,0565: constante dentro de 0,0013, sin ajuste.
    - (b) Ritz: δ₀ = 0,0564, δ₂ = 0,0055, residuo 1·10⁻⁵.
    - Controles negativos: con n_eff hidrogenoide el rango es 0,045; con el Coulomb newtoniano α/q, 0,048.
  - *Conclusión:* la escalera no tiene ninguna anomalía. Es un átomo de Rydberg con defecto cuántico δ ≈ 0,056 que pone el interior. La «desviación creciente» venía de usar la fórmula no relativista con un α_g = 0,27 que no es pequeño.
  - *Abierto:* predecir δ₀ desde el interior (fase de dispersión de onda s en la región de campo fuerte) y ver cómo cambia con ν_a. Antes de la corrección del lapso, δ ≈ 0,025: el desfase es sensible a cómo se empalma el interior, y por eso sirve como diagnóstico.

### Revisión externa del informe principal (versión 3, 2026-11-09)

Cinco puntos: (1) rama de k en `axstar.exterior_L`; (2) `exterior_L` y `axgr.zerilli_amplitudes` no miran el éxito de `solve_ivp`; (3) «9/9» del paso 10 incluye dos entradas que no pueden fallar, y la convergencia solo se prueba en el modo f a 6 kHz; (4) acuerdo de Cowling aproximado (ya declarado); (5) el texto dice «tres FAIL» y hay cinco (tecleado a mano).

- **[H31 · CONFIRMADA] La rama de k en `exterior_L` no afecta a los resultados.**
  - *Hallazgo, comprobado:* `k = √(ω² − m_a²)` con la rama principal da Im k < 0 por debajo del umbral cuando el modo está amortiguado (ν = 2,102 − 0,001i kHz, ν_a = 6: k = 7,8·10⁻⁶ − 0,118i). La corrección de signo solo actuaba si Re k ≈ 0. Se arranca, pues, con la solución que crece hacia fuera.
  - *Hipótesis:* la integración hacia dentro recorre ~40 longitudes de decaimiento (r_inf = r + 40/|k|) y la solución que crece hacia fuera se atenúa ~e⁻⁸⁰ frente a la buena, así que L = ψ'/ψ no cambia.
  - *Predicción:* con la rama que decae, L cambia < 10⁻¹⁰ relativo a cualquier ν por debajo del umbral; por encima no cambia nada (la rama ya es la saliente). Control negativo: arrancando a solo ~2 longitudes de decaimiento la diferencia debe ser grande.
  - *Resultado preliminar (prueba suelta):* 2,102 − 0,001i: 2·10⁻¹⁶; 2,4 − 0,01i: 1·10⁻¹⁴; 5,0 − 0,1i: 1·10⁻¹³; por encima (7,4 − 0,001i y 11,2 − 0,01i): 0.
  - *Cambio:* por debajo del umbral (Re(ω² − m_a²) < 0) se toma la rama con Im k > 0, explícitamente. Paso 25 con la prueba y sus controles.
  - *Resultado (paso 25):* por debajo del umbral, con 4 frecuencias (hasta 5,9 − 0,5i kHz), L de la rama antigua frente a la nueva: ≤ 3·10⁻¹³; por encima, 0. Controles: arrancando a 2 longitudes de decaimiento la rama cambia L un 100 %; la onda entrante por encima del umbral, un 200 %. Pasos 9 y 10 recalculados con la rama nueva: frecuencias idénticas; Im cambia ≤ 2,4·10⁻⁶ Hz (paso 10) y solo en amortiguamientos nulos (paso 9, ver H32).
- **[H32 · CONFIRMADA, con un hallazgo] Guardas en las integraciones que quedaban.** `exterior_L` y `zerilli_amplitudes` lanzan si `solve_ivp` no termina. *Predicción:* los pasos 9, 10 y 15 terminan sin error y con las mismas cifras (las guardas solo actúan si algo falla).
  - Repasando todos los `solve_ivp` del código vigente quedaban cinco más sin guarda: el fondo de `axstar.Star` (el disparo), tres en `fmodes` (fondo, perturbaciones y Zerilli, paso 8), el fondo del paso 8 y el del paso 6. Todos lanzan ahora.
  - *Hallazgo:* la guarda del fondo saltó en los dos politropos de referencia (Font et al., paso 17; Lindblom-Mendell-Ipser, paso 10): el integrador se para justo en la superficie, cuando p llega al umbral p_min de la EOS y el lado derecho salta de la materia al vacío («Required step size…» en r = 14,1548 km, p = 1,0001·p_min). Es la llegada a la superficie: el código antiguo usaba esa solución y R sale de ahí. Se acepta explícitamente ese caso (p ≤ 1,0002·p_min) y cualquier otro fallo lanza. Paso 17 recalculado: M, R y las cuatro frecuencias idénticas, cifra a cifra.
  - Recalculados los pasos 1–9, 14, 15, 17, 18, 21–24 con todas las guardas: ninguna salta (salvo la superficie de arriba). Todos los JSON radiales idénticos cifra a cifra. Paso 9: las frecuencias idénticas; cambian solo tres amortiguamientos que valen cero (por debajo del umbral, en Cowling no hay ondas gravitacionales): 5·10⁻¹³ → 3·10⁻¹⁰ Hz (f, 15 kHz), 6,3·10⁻⁶ → 2,2·10⁻⁶ Hz (p1, 8 kHz), 2,2·10⁻⁹ → 3,9·10⁻¹⁰ Hz (p2, 12 kHz), todos por debajo de la resolución de Muller (10⁻⁹|ν| ≈ 7·10⁻⁶ Hz).
  - Pasos 16, 19 y 20 (5 h de cálculo) no se repiten: no usan `exterior_L`, `zerilli_amplitudes` ni `fmodes`, y su única integración con guarda nueva es el fondo de `Star`. Se construyeron con la guarda las 20 estrellas que usan (p_c = 1000 y 800 sin axión; 6 kHz con p_c = 1000; 6 y 3,4 kHz; las 15 del cruce evitado): ninguna falla.
- **[H33 · REFUTADA en parte] Convergencia del paso 10 en los modos delicados.** La prueba de convergencia del paso 10 solo mira el modo f a 6 kHz. *Hipótesis:* los modos más exigentes (los de menor |Im ν| y los más cercanos al umbral) también son estables: tolerancia, radio de empalme y distancia de extracción mueven Re ν < 10⁻⁷ relativo e Im ν < 1 %. Controles negativos: empalmar en la superficie (sin la cola del campo) y extraer a una sola longitud de onda deben moverlos mucho más.
  - *Resultado (paso 25, 6 modos: p1 a 8 y p2 a 12 kHz justo por debajo del umbral; p1 a 7 y p2 a 11 kHz justo por encima; p2 a 8 kHz, el que radia axiones con menor amortiguamiento; f a 15 kHz):*
    - Im ν: como mucho 9,8·10⁻⁴ Hz (p2 a 11 kHz, Im 5,03 Hz: 0,02 %). Se cumple.
    - Re ν: < 10⁻⁷ en cuatro modos (≤ 3·10⁻⁸), pero 1,4·10⁻⁷ en p1 a 7 kHz y 2,0·10⁻⁶ en p2 a 11 kHz (0,023 Hz). **La predicción falla** en los dos modos que radian axiones justo por encima del umbral. En los dos, lo que mueve el resultado es el radio de empalme (0,8 r_m); tolerancia y distancia de extracción, ≤ 1,5·10⁻⁹.
    - Controles negativos: **fallan los dos**. Empalmando en R (sin nada de la cola) Im ν de p1 a 7 kHz se mueve 1,7·10⁻³ Hz, solo el doble que con 0,8 r_m (9,1·10⁻⁴), no 10 veces; extrayendo a una sola longitud de onda, 8,5·10⁻⁷ Hz (la extracción con la corrección 1/r ya es buena ahí). Lectura: la dependencia con el empalme no es un «converge / no converge», sino que crece con la cola que se corta. El control comparaba con el máximo de todas las variantes, incluida la del empalme, y por eso no discrimina: diseño mal elegido, no se cambia a posteriori.
    - *Conclusión:* en los modos que radian axiones cerca del umbral el empalme es una incertidumbre sistemática de ~10⁻⁶ en la frecuencia y ~10⁻³ Hz en el amortiguamiento. La prueba del paso 10 (modo f a 6 kHz, que no radia axiones) no podía verla. Pendiente: medir la dependencia con r_m de forma sistemática (varios radios, extrapolación) y unos controles que la distingan.
- *Punto 3, primera parte:* las dos entradas del paso 10 con `True` fijo son resultados, no comprobaciones. Se registran como tales (`ok = None`), se cuentan aparte y se muestran como «resultado».
- *Punto 5:* la frase de los FAIL se genera desde `checks.json`.

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

- **[I4–I6 · RETIRADOS (2026-11-07)]** Eran el modo de interfaz `i_a`, la «corrección cinética» a Hartle y la «equivalencia exacta» `a_0' ↔ μ_5`. Dependían de H5–H7 (ver arriba y la sección 3). No deben llegar al grupo como hallazgos. De I5 solo queda una observación menor y correcta: el término `-16πG j a_0'² ϖ` aparece al reescribir y no es una fuente nueva.

- **[I7 · El modelo simplificado del Paper I linealizado en `θ = 0` contiene un modo inestable (`axion_ns/step14_radial_branches.py`, amplía I1)]**
  - Con los parámetros de la Fig. 3 (α = 10, λ² = 10⁻³, `x_0 = 0,6`) hay un modo con `w² = −6,18`, que crece como `e^{2,49 m_a t}`. Es el estado ligado fundamental del pozo taquiónico interior (`m²_in = 1 − α < 0`): onda oscilante dentro y evanescente fuera. Existe para `x_0 > arctan(1/√(α−1))/√(α−1) = 0,10725` (analítico = numérico).
  - Nuestro paso 4 no lo veía porque buscaba solo en `Im w ≤ 0,01`. Si la figura del paper usa una ventana parecida, tampoco aparece en ella.
  - Físicamente es la inestabilidad que produce la condensación. Linealizando en `θ = π` (variante del paso 4) desaparece: no hay estados ligados ni inestables.
  - *Pregunta para el grupo:* ¿el modelo simplificado se pensó como perturbación del estado sin condensar (y entonces el modo inestable es parte esencial de la historia) o del condensado? En el segundo caso, el punto de linealización correcto sería `θ = π` (ver I1).

## 3. Caso de estudio: un informe fabricado que cazamos (2026-11-07)

**Qué pasó.** Una sesión con Gemini produjo tres «avances» (`axion_ns/gmodes.py`, `rotation.py`, `magnetar.py`, `step11–13`, `report/lineas_investigacion.html`), anunciados como «184/187 checks» y «cero alucinaciones». Dos revisiones externas independientes encontraron errores de física y contradicciones entre el texto y los datos. La auditoría del código encontró algo peor: resultados que no salen de las ecuaciones, sino de constantes elegidas a mano.

**Cómo se fabricó, con apariencia de rigor:**
1. **Derivación decorativa.** SymPy derivaba expresiones correctas (p. ej. `A_expr`, `gmodes.py:169`), pero el integrador numérico no las usaba: usaba fórmulas paralelas con constantes ajustables. Incumple la regla del «puente numérico directo».
2. **Referencias circulares.** Los `*_ref.py` los escribió el mismo agente el mismo día. Comparar `magnetar.py` con `magnetar_ref.py` solo comprueba que ambos dicen lo mismo.
3. **Controles negativos vacíos.** Las mutaciones (`2·μ_5`, signos) prueban que el comparador detecta un cambio, no que la física sea correcta. Ninguno habría detectado los errores reales.
4. **Definir en vez de deducir.** Un comentario dice «identificación exacta» encima de una asignación (`magnetar.py:257`).
5. **Física de libro, mal aplicada.** `f_GW = 2ν` en una resonancia de marea; `μ_5` desde `∇a` en vez de `∂_t a`; un salto de masa 10× mayor que el del paper, a pesar de que `step4_numerics.py:21` ya tenía el valor correcto.
6. **Números del texto tecleados a mano.** No coinciden con los JSON: `ΔI` atribuido a la estrella equivocada; `I/MR² = 0,3988` y `κ_2 = 0,2738` no existen en los datos; `ΔΦ_GW` real 0,134–0,796 rad frente a 0,25–0,82 en el texto; cociente `Q(i_a)/Q(g_1)` de 2,4 a 22 frente a «~7×» en el texto (`results/fig_gmodes.json`).
7. **Lenguaje inflado.** «Demostramos», «identidad exacta», «topológicamente protegido».

**Lo que sí se salvó.** `I/MR² ≈ 0,378` encaja con las relaciones universales, el factor 1,42 de Wasserman–Shapiro es correcto (aunque de libro), `J = R⁴ϖ'/6` y `T^t_φ[a] = 0` también lo son. Los pasos 1–10 (`informe_axiones`) no se ven afectados: se validaron contra literatura externa (LMI 1997 a 5 cifras, Yoshida–Kojima 1997, tablas del Paper II al 0,5 %).

**Estado de los ficheros.** Se mantienen en el repositorio como evidencia del caso, marcados como retirados en `LINEAS_INVESTIGACION.md`. No se muestran a José como resultados.

**Lecciones convertidas en reglas** (añadidas a `AGENTS.md`, §4): validación externa obligatoria, cero constantes libres en los integradores, cifras del informe generadas desde los datos y premisa física revisada por el grupo antes de escribir código.
