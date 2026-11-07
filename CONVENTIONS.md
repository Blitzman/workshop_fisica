# Convenciones físicas y de notación del grupo (Gómez-Bañón, Pnigouras & Pons)

## 1. Métrica, signatura y tensores
- **Signatura:** `(+1, -1, -1, -1)`.
- **Métrica estática de fondo:** `ds² = e^{2Φ(r)} dt² - e^{2λ(r)} dr² - r² (dθ² + sin²θ dφ²)`, con `e^{-2λ} = 1 - 2GM(r)/r`.
- **Convenio de curvatura (`axion_ns/geometry.py`):**
  - `Γ^a_{bc} = ½ g^{ad}(∂_b g_{dc} + ∂_c g_{db} - ∂_d g_{bc})`
  - `R^a_{bcd} = ∂_c Γ^a_{db} - ∂_d Γ^a_{cb} + Γ^a_{ce} Γ^e_{db} - Γ^a_{de} Γ^e_{cb}`, `R_{bd} = R^a_{bad}`
  - Ecuaciones de Einstein: `G_{μν} = 8πG T_{μν}`
- **Tensor energía-momento del fluido:** `T_{μν} = (ε + p) u_μ u_ν - p g_{μν}`.

## 2. Perturbaciones no radiales (Cowling y Regge-Wheeler / Lindblom-Detweiler)
- **Dependencia temporal y armónica:** `e^{-iωt} Y_{lm}(θ, φ)` (modo polar axisimétrico `m = 0`; las ecuaciones radiales no dependen de `m`).
- **Factor radial extraído `(r/R)^l` (`axion_ns/nonradial.py`):**
  - Axión: `a = a_0(r) + f_a (r/R)^l δθ(r) Y_{lm} e^{-iωt}`
  - Bariones: `n_b = n_0(r) + (r/R)^l δn_b(r) Y_{lm} e^{-iωt}`
  - Desplazamiento fluido: `ξ^r = (r/R)^l e^{-λ} (W(r)/r) Y_{lm} e^{-iωt}`, `ξ^θ = -(r/R)^l (V(r)/r²) ∂_θ Y_{lm} e^{-iωt}` (¡ojo al signo menos en `ξ^θ`!).
- **Métrica perturbada (convenio Lindblom-Detweiler 1983, `axion_ns/fullgr.py`):**
  `ds² = e^{2Φ}(1 + H_0 rY) dt² - 2iωr H_1 rY dt dr - e^{2λ}(1 - H_2 rY) dr² - r²(1 - K rY) dΩ²`, donde `rY = (r/R)^l Y_{lm} e^{-iωt}`.
- **Variable de presión lagrangiana:** `X = -e^Φ Δp = -e^Φ (δp + ξ^r p')`.
- **Modos cuasinormales:** `ω = ω_R + i ω_I` con `ω_I < 0` para modos amortiguados (`e^{-iωt} = e^{-|ω_I|t} e^{-iω_R t}`, `1/τ = 2π |Im ν|`).

## 3. Ecuación de estado (EOS) y unidades
- **Unidades simbólicas:** `c = ħ = 1`, `G` explícita.
- **Unidades en integradores numéricos (`axion_ns/axstar.py`, `axion_ns/axgr.py`):** Geométricas en kilómetros (`G = c = 1`, densidades y presiones en `km⁻²` vía `KGEO = 1.32379e-6 km⁻² / (MeV fm⁻³)`).
- **Politropos de calibración (LMI 1997 / Yoshida-Kojima 1997):** En `p = κ ρ²`, `ρ` denota **densidad de energía `ε`**, NO densidad de masa en reposo `m_b n_b`.
- **BSk26 y axión:** `ε(n, a) = ε_0(n) + n [m*(a) - m_N]`, por lo que `p = n ∂ε/∂n - ε = p_0(n)` no depende explícitamente de `a`. Escala axiónica fijada por `E_0 = ε m_π² f_π² = n_c σ_N`.
