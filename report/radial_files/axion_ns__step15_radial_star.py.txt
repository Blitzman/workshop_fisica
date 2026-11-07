"""Paso 15: modos radiales de la estrella BSk26 con condensado en relatividad general completa (hipótesis H11).

Ecuaciones: axion_ns/radial.py (derivadas desde la acción, sin transcribir). Las semillas de los buscadores de raíces
son orientativas (exploración previa o fórmulas): cada raíz se acepta solo si el determinante se anula y su carácter
(fracción de energía en el fluido) es el esperado.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from . import axgr as AG
from . import axstar as AX
from . import radial as RD
from .pipeline import report

OUT = Path(__file__).resolve().parent.parent / 'results'
MSUN_KM = 1.476625


def exterior_vs_37(drop_mass=False):
    """Ecuación exterior derivada con theta_0 = 0 frente a la Ec. (3.7) del Paper I (verificada en el paso 2)."""
    Ae = RD.exterior_matrix()
    ma = AX.omega_of(6.0)
    E0 = AX.E0 * AX.KGEO
    fa = np.sqrt(AX.BETA * E0 / (4 * ma**2))
    M, r = 2.17, 20.0
    f = 1 - 2 * M / r
    w = AX.omega_of(4.2)
    P1 = M / (r * (r - 2 * M))
    P2 = -2 * (r - M) * M / (r * (r - 2 * M))**2
    G2 = 0.0 if drop_mass else AX.g2(0.0)
    bg = (0.5 * np.log(f), P1, P2, -0.5 * np.log(f), -P1, 0, 0, 0.0, 0.0, 0, 0, 1.0, fa, 0.0, 0.0, E0 * G2, 0.0, -E0 * G2)
    A = np.asarray(Ae(r, w, *bg), dtype=complex)
    want = np.array([(ma**2 - w**2 / f) / f, -(2 / r) * (1 - M / r) / f])
    return float(max(abs(A[1, 0] - want[0]) / abs(want[0]), abs(A[1, 1] - want[1]) / abs(want[1])))


# ------------------------------------------------------------------ (a) sin axión: punto de retorno
def wF2(st, lo=-3e-3, hi=3e-3):
    g = lambda w2: RD.interior(st, np.sqrt(complex(w2)), 'fluid').y[1, -1].real
    while g(lo) * g(hi) > 0:
        lo, hi = lo * 2, hi * 2
    return brentq(g, lo, hi, xtol=1e-12)


def turning_point():
    M = lambda p0: -AX.Star(None, p0=p0).M
    pm = minimize_scalar(M, bounds=(700, 1100), method='bounded', options={'xatol': 0.5}).x
    pz = brentq(lambda p0: wF2(AX.Star(None, p0=p0)), 800, 1000, xtol=0.5)
    return float(pm), float(pz), float(-M(pm) / MSUN_KM)


# ------------------------------------------------------------------ (b, c) con axión
def real_root(st, lo, hi, rm):
    return brentq(lambda nu: RD.det(st, AX.omega_of(nu), rm).real, lo, hi, xtol=1e-9)


def complex_root(st, seed, rm):
    f = lambda nu: RD.det(st, AX.omega_of(nu), rm)
    z = RD.muller(f, seed)
    return z, abs(f(z))


def run():
    res = []
    out = {}
    # 1. Ecuación exterior = Ec. (3.7) y control negativo con el error real que cazamos (sin término de masa)
    e37, e37_bad = exterior_vs_37(), exterior_vs_37(drop_mass=True)
    res.append((f'Ecuación exterior derivada (theta_0 = 0) = Ec. (3.7) del Paper I: error relativo {e37:.1e}', e37 < 1e-12, ''))
    res.append((f'Control negativo 1 (sin término de masa, el error real del 2026-11-07): error {e37_bad:.2f}', e37_bad > 0.1, ''))

    # 2. Sin axión: omega_F^2 = 0 en la masa máxima (punto de retorno)
    pm, pz, Mmax = turning_point()
    out['turning'] = {'p0_Mmax': pm, 'p0_wF2_zero': pz, 'Mmax_msun': Mmax}
    res.append((f'Sin axión (BSk26): masa máxima {Mmax:.4f} M☉ en p_c = {pm:.1f} MeV/fm³; ω_F² = 0 en p_c = {pz:.1f} '
                f'(diferencia {abs(pm - pz) / pm * 100:.2f} %)', abs(pm - pz) / pm < 0.01, 'criterio de Chandrasekhar/punto de retorno'))

    # 3. Exterior analítico en el rayo complejo frente a la recta real (modo muy amortiguado)
    st6 = AX.Star(6.0)
    rm6 = AG.r_match(st6)
    Mt, fa_v = st6.sol.y[0, -1], np.sqrt(st6.fa2)
    z0, h = AX.omega_of(6.585 - 0.563j), AX.omega_of(1e-3)
    Lc = lambda om: (lambda p: p[1] / p[0])(RD.exterior_psi(om, Mt, rm6, st6.ma, fa_v, st6.E0))
    Lr = lambda om: AX.exterior_L(st6, 0, om, r_stop=rm6)
    cr = lambda L: abs((L(z0 + h) - L(z0)) / h - (L(z0 + 1j * h) - L(z0)) / (1j * h)) / abs((L(z0 + h) - L(z0)) / h)
    crc = cr(Lc)
    rel = lambda z: abs(Lr(AX.omega_of(z)) - Lc(AX.omega_of(z))) / abs(Lc(AX.omega_of(z)))
    d_strong, d_weak, d_below = rel(6.585 - 0.563j), rel(8.132 - 0.0036j), rel(5.5)
    crr = d_strong
    res.append((f'Exterior en el rayo complejo: cumple Cauchy-Riemann (desajuste {crc:.1e}). El de la recta real (axstar.exterior_L) '
                f'se desvía {d_strong * 100:.0f} % en la resonancia muy amortiguada, {d_weak:.0e} en un modo poco amortiguado y '
                f'{d_below:.0e} bajo el umbral', crc < 1e-2 and d_strong > 0.5 and d_weak < 1e-3 and d_below < 1e-8,
                'la onda saliente crece con r y la integración hacia dentro es inestable; los pasos 9-10 (Im de pocos Hz) no se ven afectados'))

    # 4. nu_a = 6 kHz: modo F ligado (frecuencia real) y de carácter fluido
    F6 = real_root(st6, 3.40, 3.50, rm6)
    fF6 = RD.fluid_fraction(st6, AX.omega_of(F6), rm6)
    imrel = abs(RD.det(st6, AX.omega_of(4.5), rm6).imag) / abs(RD.det(st6, AX.omega_of(4.5), rm6))
    res.append((f'ν_a = 6 kHz: modo F en {F6:.4f} kHz < ν_a, frecuencia real (|Im det|/|det| = {imrel:.0e} en el eje), '
                f'fracción fluida {fF6:.3f}', imrel < 1e-10 and fF6 > 0.9, 'H11(b)'))

    # 5. Escalera hidrogenoide: n_eff entero
    alpha_g = Mt * st6.ma
    ladder = []
    for n in range(2, 7):
        pred = 6.0 * (1 - alpha_g**2 / (2 * n * n))
        dn = 6.0 * alpha_g**2 / n**3 * 0.3
        z = real_root(st6, pred - dn, min(pred + dn, 6.0 * (1 - 1e-6)), rm6)
        ladder.append((n, z, alpha_g / np.sqrt(2 * (1 - z / 6.0)), RD.fluid_fraction(st6, AX.omega_of(z), rm6) if n <= 3 else None))
    ok_l = all(abs(ne - n) < 0.05 for n, _, ne, _ in ladder) and all(fr < 0.05 for *_, fr in ladder if fr is not None)
    res.append((f'Átomo gravitatorio (α_g = {alpha_g:.4f}): n_eff = ' + ', '.join(f'{ne:.3f}' for _, _, ne, _ in ladder)
                + ' para n = 2…6; carácter axiónico', ok_l, 'H11(c)'))

    # 6. Por encima del umbral: H1 (fluido, poco amortiguado) y resonancia axiónica (muy amortiguada)
    H1, rH1 = complex_root(st6, 8.13 - 0.004j, rm6)
    fH1 = RD.fluid_fraction(st6, AX.omega_of(H1), rm6)
    Ra, rRa = complex_root(st6, 6.6 - 0.55j, rm6)
    fRa = RD.fluid_fraction(st6, AX.omega_of(Ra), rm6)
    res.append((f'ν_a = 6 kHz: H1 = {H1.real:.4f} {H1.imag:+.4f}i kHz (fracción fluida {fH1:.3f}); resonancia axiónica '
                f'{Ra.real:.3f} {Ra.imag:+.3f}i kHz (fracción fluida {fRa:.3f})',
                rH1 < 1e-11 and rRa < 1e-11 and H1.imag < 0 and fH1 > 0.9 and Ra.imag < -0.3 and fRa < 0.05, 'H11(b) y H10'))
    d_ok = abs(RD.det(st6, AX.omega_of(H1), rm6, sheet=-1)) / max(abs(RD.det(st6, AX.omega_of(H1 + 0.05), rm6)), 1e-300)
    res.append((f'Control negativo 2 (onda entrante): en la raíz H1 el determinante deja de anularse (|det|/escala = {d_ok:.1e})',
                d_ok > 1e-3, ''))

    # 7. La transición: cruce evitado F – átomo n = 1 y paso de F por el umbral
    st35, st34 = AX.Star(3.5), AX.Star(3.4)
    rm35, rm34 = AG.r_match(st35), AG.r_match(st34)
    b35 = [real_root(st35, 3.38, 3.43, rm35), real_root(st35, 3.44, 3.48, rm35)]
    f35 = [RD.fluid_fraction(st35, AX.omega_of(z), rm35) for z in b35]
    F34, rF34 = complex_root(st34, 3.42 - 0.001j, rm34)
    fF34 = RD.fluid_fraction(st34, AX.omega_of(F34), rm34)
    res.append((f'ν_a = 3,5 kHz: dos ligados {b35[0]:.4f} y {b35[1]:.4f} kHz con carácter mezclado '
                f'({f35[0]:.2f} y {f35[1]:.2f} de fluido): cruce evitado F – átomo n = 1',
                all(0.1 < x < 0.9 for x in f35), ''))
    res.append((f'ν_a = 3,4 kHz: F ya por encima del umbral, Re ν = {F34.real:.4f} kHz e Im ν = {F34.imag * 1e3:.2f} Hz '
                f'(fracción fluida {fF34:.3f})',
                rF34 < 1e-11 and F34.real > 3.4 and F34.imag < 0 and fF34 > 0.9, 'el modo fluido pasa a cuasi-normal al cruzar ν_a'))

    # 8. Revisión externa (2026-11-08): estabilidad radial (H12) y robustez numérica de los amortiguamientos (H13)
    stab = {}
    for lab, st_, rm_ in (('6', st6, rm6), ('3.4', st34, rm34)):
        gs = np.geomspace(0.01, 10, 40)
        v = [RD.det(st_, AX.omega_of(1j * g), rm_) for g in gs]
        stab[lab] = {'sign_changes': int(sum(v[i].real * v[i + 1].real < 0 for i in range(len(gs) - 1))),
                     'im_rel': float(max(abs(x.imag) / abs(x) for x in v))}
    res.append(('Estabilidad radial: sin raíces en ω = iγ, γ ∈ [0,01; 10] kHz, para ν_a = 6 y 3,4 kHz (' +
                ', '.join(f'{k}: {v_["sign_changes"]} cambios de signo' for k, v_ in stab.items()) + ')',
                all(v_['sign_changes'] == 0 and v_['im_rel'] < 1e-10 for v_ in stab.values()),
                'H12; un barrido de signo puede no ver pares de raíces muy próximas'))
    robust = {}
    for lab, st_, rm_, z0 in (('H1 (6 kHz)', st6, rm6, H1), ('F (3,4 kHz)', st34, rm34, F34)):
        rr = np.linspace(st_.R, rm_, 4000)
        th = np.abs(st_.sol.sol(rr)[3])
        rows = {'base': z0}
        for name, kw in (('rtol 1e-12', {'rtol': 1e-12}), ('r0 1e-4 km', {'r0': 1e-4}), ('exterior ×1,5', {'n_decay': 60.0})):
            rows[name] = RD.muller(lambda nu: RD.det(st_, AX.omega_of(nu), rm_, **kw), z0)
        for tol in (1e-5, 1e-6):
            r2 = float(rr[np.argmax(th < tol)])
            rows[f'empalme |θ0|<{tol:g} ({r2:.0f} km)'] = RD.muller(lambda nu: RD.det(st_, AX.omega_of(nu), r2), z0)
        robust[lab] = {k: [v_.real, v_.imag] for k, v_ in rows.items()}
    num = max(abs(v_[1] - r_['base'][1]) / abs(r_['base'][1]) for r_ in robust.values()
              for k, v_ in r_.items() if k in ('rtol 1e-12', 'r0 1e-4 km', 'exterior ×1,5'))
    emp = {lab: abs(list(r_.values())[-1][1] - r_['base'][1]) / abs(r_['base'][1]) for lab, r_ in robust.items()}
    res.append((f'Robustez de Im ν: rtol, r0 y longitud exterior cambian como mucho un {num * 100:.3f} %; el radio de empalme '
                '(último paso hacia fuera) un ' + ', '.join(f'{v_ * 100:.2f} % en {k}' for k, v_ in emp.items()),
                num < 5e-4 and all(v_ < 0.01 for v_ in emp.values()), 'H13; el radio de empalme es la incertidumbre dominante'))

    out.update({
        'stability': stab, 'robustness': robust,
        'nu6': {'M_msun': st6.M / MSUN_KM, 'R_km': st6.R, 'alpha_g': alpha_g, 'F': F6, 'F_fluid': fF6,
                'ladder': [{'n': n, 'nu': z, 'n_eff': ne} for n, z, ne, _ in ladder],
                'H1': [H1.real, H1.imag], 'H1_fluid': fH1, 'axion_res': [Ra.real, Ra.imag], 'axion_res_fluid': fRa},
        'nu35': {'bound': b35, 'fluid': f35}, 'nu34': {'F': [F34.real, F34.imag], 'fluid': fF34},
        'exterior': {'cr_complex_ray': crc, 'cr_real_line': crr, 'old_vs_new_strong': d_strong, 'old_vs_new_weak': d_weak},
        'checks': [{'name': n_, 'ok': bool(o_), 'note': t_} for n_, o_, t_ in res],
    })
    OUT.mkdir(exist_ok=True)
    (OUT / 'radial_star.json').write_text(json.dumps(out, indent=1))
    return report('Paso 15: modos radiales de la estrella con condensado (GR completa)', res)


if __name__ == '__main__':
    run()
