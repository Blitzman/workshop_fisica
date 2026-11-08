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
  - **[H17 · CONFIRMADA] ARTEFACTO RETIRADO: la «rama axiónica» de la estrella realista no existe.** Eran las resonancias muy amortiguadas y casi equiespaciadas (6,585 − 0,563 i; 7,304 − 0,748 i; 8,142 − 0,884 i; 9,062 − 0,985 i; 10,039 − 1,058 i kHz).
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
    - Sigue pendiente y la reclaman las tres revisiones: el principio del argumento en el semiplano superior para la estrella realista. También una comprobación independiente de las ecuaciones radiales en GR (por ejemplo, frente a un código publicado de modos radiales sin axión, más allá del punto de retorno).
    - El defecto cuántico negativo de la escalera (n_eff = 5,967 en n = 6) sigue sin explicar.
    - En `build_radial.py` quedan cifras tecleadas en el recuadro de la corrección (~0,16 %, ~0,5 %, 1,6·10⁻³). Para el informe radial faltan además: la estabilidad con el principio del argumento en el semiplano superior (crítica A), el barrido denso en el umbral (D), redactar como «no detectamos inestabilidades» y no como «es estable», y separar mejor lo demostrado de lo que es hipótesis.
  - `axstar.Phi_shift` ajusta Φ a Schwarzschild en R con `m(R)`, pero la cola del condensado fuera de la estrella lleva un 0,5 % de la masa. Lejos, `Φ − ½ ln(1 − 2M_t/r) = +1,6·10⁻³` (a 6 y a 3,4 kHz), en vez de 0.
  - *Predicción:* normalizando con `M_t` al final de la integración del fondo, (i) la dependencia de `Im ν` con el radio de empalme cae a ≲0,1 % y (ii) todas las frecuencias con axión se desplazan ~0,16 %.
  - *Afecta a código compartido:* los pasos 9 (Cowling, tablas del Paper II) y 10 (GR no radial) usan el mismo `Phi_shift`.

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
