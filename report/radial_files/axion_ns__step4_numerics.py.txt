"""Paso 4: comprobaciones numéricas de las cifras del paper y reproducción de la Fig. del modelo simplificado."""
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import constants as C

from .pipeline import report

OUT = Path(__file__).resolve().parent.parent / 'results'
HBARC = C.hbar * C.c / C.e * 1e9  # MeV fm  (197.327)


def constants_checks():
    res = []
    z = 0.48
    beta = 4 * z / (1 + z)**2
    res.append((f'beta = 4z/(1+z)^2 = {beta:.3f} (paper: 0.88)', abs(beta - 0.88) < 0.01, ''))
    sigma = 50.0
    dm = sigma * (1 - np.sqrt(1 - beta))
    res.append((f'Reducción de masa sigma_N [1-g(pi f_a)] = {dm:.1f} MeV (paper: ~32 MeV)', abs(dm - 32) < 1.5, ''))

    nu, R = 3e3, 12e3  # 3 kHz, 12 km
    E_peV = C.h * nu / C.e * 1e12
    res.append((f'm_a = 3 kHz <-> {E_peV:.1f} peV (paper: 12 peV)', abs(E_peV - 12) < 0.5, 'E = h nu'))
    maR_paper = nu * R / C.c
    maR_true = 2 * np.pi * nu * R / C.c
    res.append((f'm_a R para 3 kHz y R=12 km = {maR_true:.2f} (paper: 0.12)', abs(maR_true - 0.12) < 0.02,
                f'el paper usa nu R/c = {maR_paper:.2f}; m_a c^2/hbar = 2 pi nu  -> factor 2 pi'))

    rows = []
    for mpi, fpi in [(134.98, 92.1), (139.57, 92.1), (139.57, 93.0)]:
        m2f2 = (mpi * fpi)**2 / HBARC**3          # m_pi^2 f_pi^2 en MeV/fm^3
        rows.append((mpi, fpi, m2f2, beta * m2f2 / 4, m2f2 / sigma))
    lo, hi = min(r_[3] for r_ in rows), max(r_[3] for r_ in rows)
    res.append((f'm_a^2 f_a^2 = eps beta m_pi^2 f_pi^2 / 4 = eps x ({lo:.1f}-{hi:.1f}) MeV/fm^3 (paper: eps x 13)',
                lo <= 13 <= hi, f'm_pi^2 f_pi^2 = {rows[0][2]:.1f}-{rows[-1][2]:.1f} MeV/fm^3; '
                'el paper siguiente (2610.37308) da eps m_pi^2 f_pi^2 ~ 2 MeV/fm^3 para eps=0.1, coherente con esto'))
    lam2 = np.array([r_[3] for r_ in rows]) / (500 * 0.3)
    res.append((f'lambda^2 con (eps+p)=500 MeV/fm^3, c_s^2=0.3: eps x ({lam2.min():.3f}-{lam2.max():.3f}) (paper: eps x 0.02)',
                lam2.min() <= 0.02 <= lam2.max(),
                f'incluso con el 13 del paper: 13/150 = {13/150:.3f}; las dos cifras del paper no son coherentes entre sí'))
    nc = rows[0][4]
    res.append((f'n_c(eps=1) = m_pi^2 f_pi^2/sigma_N = {nc:.2f} fm^-3 = {nc/0.16:.1f} n_sat (paper: ~2 n_sat)',
                1.5 < nc / 0.16 < 3.5, ''))
    return res, beta


# ------------------------------------------------------------------ modelo simplificado
def f_toy(w, lam2, alpha, cs2, x0, mass2_in=None, entire=True):
    """Ec. (4.13). mass2_in = masa efectiva^2 interior del axión en unidades de m_a^2 (paper: 1 - alpha).

    Con entire=True devuelve f cos(k+ x0) cos(k- x0), que es entera: en la forma con tangentes cada modo
    fluido es un cero pegado a un polo de tan(k+ x0) (separación O(lambda^2)), y los buscadores de raíces
    y los mapas de |f| no lo ven. Parece ser también la forma que se representa en la figura del paper.
    """
    if mass2_in is None:
        mass2_in = 1 - alpha
    ka2, ks2 = w**2 - mass2_in, w**2 / cs2
    A = alpha**2 * lam2
    D = (ka2 - ks2 + A)**2 + 4 * A * ks2
    sD = mp.sqrt(D)
    kp2, km2 = (ka2 + ks2 + A + sD) / 2, (ka2 + ks2 + A - sD) / 2
    kp, km = mp.sqrt(kp2), mp.sqrt(km2)
    kout = w * mp.sqrt(1 - 1 / w**2)  # rama ~ w para |w| grande: onda saliente con e^{-i w t}
    if not entire:
        return (ka2 - km2) * kp * mp.tan(kp * x0) - (ka2 - kp2) * km * mp.tan(km * x0) + 1j * kout * sD
    cp, cm = mp.cos(kp * x0), mp.cos(km * x0)
    return (ka2 - km2) * kp * mp.sin(kp * x0) * cm - (ka2 - kp2) * km * mp.sin(km * x0) * cp + 1j * kout * sD * cp * cm


def landscape(pars, re=(0.05, 30), im=(-6, 0.25), N=(600, 250), **kw):
    X, Y = np.linspace(*re, N[0]), np.linspace(*im, N[1])
    F = np.empty((N[1], N[0]))
    for j, y in enumerate(Y):
        for i, x in enumerate(X):
            F[j, i] = float(mp.log10(abs(f_toy(mp.mpc(x, y), *pars, **kw)) + 1e-300))
    return X, Y, F


def find_roots(pars, X, Y, F, **kw):
    """Semillas = mínimos locales de log|f| en la malla; refinado con Muller (mpmath)."""
    seeds = []
    for j in range(1, F.shape[0] - 1):
        for i in range(1, F.shape[1] - 1):
            if F[j, i] == F[j - 1:j + 2, i - 1:i + 2].min():
                seeds.append(complex(X[i], Y[j]))
    roots = []
    for s in seeds:
        try:
            z = mp.findroot(lambda w: f_toy(w, *pars, **kw), (s, s + 0.01, s + 0.01j), solver='muller', tol=1e-20)
        except (ValueError, ZeroDivisionError):
            continue
        z = complex(z)
        if abs(f_toy(mp.mpc(z), *pars, **kw)) < 1e-8 and X[0] <= z.real <= X[-1] and Y[0] <= z.imag <= 0.01 \
                and all(abs(z - q) > 1e-4 for q in roots):
            roots.append(z)
    return sorted(roots, key=lambda z: z.real)


def plot(panels, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 4.4), sharey=True, constrained_layout=True)
    for ax, (title, X, Y, F, roots) in zip(np.atleast_1d(axes), panels):
        im = ax.pcolormesh(X, Y, F, cmap='Blues_r', vmin=-1.5, vmax=6, shading='auto', rasterized=True)
        if roots:
            ax.plot([z.real for z in roots], [z.imag for z in roots], 'o', ms=8, mfc='none', mec='#d9480f', mew=2,
                    label='raíces (Muller)')
            ax.legend(loc='lower right', fontsize=8, frameon=False)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel(r'Re $\bar\omega$')
    np.atleast_1d(axes)[0].set_ylabel(r'Im $\bar\omega$')
    fig.colorbar(im, ax=axes, label=r'$\log_{10}|f\,\cos(k_+\bar x_0)\cos(k_-\bar x_0)|$')
    fig.savefig(path, dpi=150)
    from .figdata import save
    save('toy', [{'x': [float(X[0]), float(X[-1])], 'y': [float(Y[0]), float(Y[-1])],
                  'nx': len(X[::2]), 'ny': len(Y[::2]),
                  'f': [int(round(10 * v)) for v in np.clip(F[::2, ::2], -1.5, 6).ravel()],
                  'roots': [[round(z.real, 4), round(z.imag, 4)] for z in roots]}
                 for _, X, Y, F, roots in panels])


def run():
    res, beta = constants_checks()
    OUT.mkdir(exist_ok=True)
    pars = (1e-3, 10.0, 0.3, 0.6)  # lambda^2, alpha, c_s^2, x0 de la figura del paper
    X, Y, F = landscape(pars)
    roots = find_roots(pars, X, Y, F)
    fluid = [z for z in roots if z.imag > -0.5]
    ax_ = [z for z in roots if z.imag <= -0.5]
    # modos fluidos ~ desacoplados: cos(k_s x0) = 0  ->  w = c_s (n+1/2) pi / x0
    pred = [np.sqrt(0.3) * (k + 0.5) * np.pi / 0.6 for k in range(len(fluid))]
    ok = len(fluid) >= 8 and all(abs(z.real - p_) < 0.3 for z, p_ in zip(fluid, pred))
    res.append((f'Fig. 3 reproducida: {len(fluid)} modos fluidos (Im~0) + {len(ax_)} modos axiónicos', ok and len(ax_) >= 4,
                'axiónicos: ' + ', '.join(f'{z.real:.2f}{z.imag:+.2f}i' for z in ax_[:5])))

    # variante consistente con theta = pi: sin acoplo lineal (dm*/da = 0) y masa efectiva positiva
    m2pi = (10.0 - 1) / np.sqrt(1 - beta)
    pars_pi = (0.0, 10.0, 0.3, 0.6)
    Xp, Yp, Fp = landscape(pars_pi, mass2_in=m2pi)
    roots_pi = find_roots(pars_pi, Xp, Yp, Fp, mass2_in=m2pi)
    fl_pi = [z for z in roots_pi if abs(z.imag) < 1e-8]
    ax_pi = [z for z in roots_pi if abs(z.imag) >= 1e-8]
    res.append((f'Variante theta=pi (m_eff^2 = {m2pi:.1f} m_a^2, sin acoplo): {len(fl_pi)} modos fluidos SIN '
                f'amortiguar + {len(ax_pi)} axiónicos', True, 'axiónicos: ' + ', '.join(f'{z.real:.2f}{z.imag:+.2f}i' for z in ax_pi[:5])))
    plot([('Modelo del paper (Ec. 4.13)\n' r'$\lambda^2=10^{-3}$, $\alpha=10$, $c_s^2=0.3$, $\bar x_0=0.6$', X, Y, F, roots),
          (r'Linealizando en $\theta=\pi$: $m_\mathrm{eff}^2=(\alpha-1)/\sqrt{1-\beta}$, $\partial_a m^*=0$', Xp, Yp, Fp,
           roots_pi)], OUT / 'toy_qnm_landscape.png')
    res.append(('Figura guardada en results/toy_qnm_landscape.png', True, ''))
    return report('Paso 4: números y figura', res)


if __name__ == '__main__':
    run()
