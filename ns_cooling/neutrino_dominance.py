"""Shade the leading neutrino process, or photons when L_gamma > total L_nu."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def save_dominance_plot(entries, path, minimum_luminosity=1e31):
    names = list(dict.fromkeys(name for _, run in entries for name, _ in run.neutrino_processes))
    names.append('none')
    colors = {name: plt.get_cmap('tab20')((2 * i) % 20 + (i // 10) % 2) for i, name in enumerate(names)}
    colors['photons'] = '#e34a33'
    fig, axes = plt.subplots(len(entries), 1, figsize=(10, 4 * len(entries)), squeeze=False)
    for ax, (label, run) in zip(axes[:, 0], entries):
        ax.set_xscale('log'); ax.set_yscale('log')
        times = run.step_times_years
        dominant = run.dominant_neutrino_process
        # Midpoints bracket discrete classifications; transitions are only resolved
        # to the accepted time-step spacing. Omit t=0 from the logarithmic plot.
        t = times[1:]
        d = np.asarray(dominant[1:], dtype=object).copy()
        # Compare photons with ALL neutrinos, at accepted sample times.
        # Do not extrapolate photon luminosity before the first output sample.
        if len(run.times_years):
            known = t >= run.times_years[0]
            photon = np.interp(t, run.times_years, run.luminosity_erg_s)
            total_neutrino = run.neutrino_luminosities[1:].sum(axis=1)
            d[known & (photon > total_neutrino)] = 'photons' 
        seen = []
        if len(t) > 1:
            edges = np.r_[t[0], np.sqrt(t[:-1] * t[1:]), t[-1]]
            boundaries = np.r_[0, np.flatnonzero(d[1:] != d[:-1]) + 1, len(d)]
            for a, b in zip(boundaries[:-1], boundaries[1:]):
                name = d[a]
                ax.axvspan(edges[a], edges[b], color=colors[name], alpha=.22, lw=0)
                if name not in seen:
                    seen.append(name)
        ax.plot(run.times_years, run.luminosity_erg_s, color='black', lw=1.8)
        ax.axhline(minimum_luminosity, color='gray', ls=':', lw=.8)
        ax.set(xlabel='Time [yr]', ylabel=r'$L_\gamma^\infty$ [erg s$^{-1}$]', title=label)
        ax.legend(handles=[Patch(facecolor=colors[n], alpha=.35, label=n.replace('_',' ')) for n in seen], fontsize=8, loc='best')
        ax.grid(alpha=.15)
    fig.suptitle('Shading: leading neutrino process; red when photons dominate', fontsize=12)
    fig.tight_layout(rect=(0,0,1,.97))
    fig.savefig(path, dpi=170)
    return fig
