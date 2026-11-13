"""Cooling curves of one star with the standard and with the updated neutrino emissivities.

    python compare_emissivities.py                       # 1.4 Msun star of input_bsk26, to 2e5 yr
    python compare_emissivities.py --input input_bsk22 --out results_emissivities

Runs the star with the standard rates (``emissivity.py``), with the updated rates of
Bottaro, Caputo and Fiorillo (``emissivity_updated.py``) and with variants that isolate
each ingredient. Writes ``emissivity_comparison.csv`` (time and luminosity of every variant),
``emissivity_comparison.png`` and prints the luminosity ratios at a set of ages.
"""

import argparse
import copy
import re
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path

import numpy as np

import emissivity as standard
import emissivity_updated as updated
from config import CoolingConfig
from cooling import run_cooling
from emissivity import NEUTRINO_PROCESSES
from input_model import load_backgrounds

ROOT = Path(__file__).resolve().parent

VARIANTS = {
    "standard": "standard rates",
    "updated": "updated rates (all corrections, electron and muon channels)",
    "updated, no muon channels": "updated rates without the muon channels of the modified Urca",
    "updated modified Urca only": "updated modified Urca (with muons), standard bremsstrahlung",
    "updated bremsstrahlung only": "updated nn, np, pp bremsstrahlung, standard modified Urca",
    "standard, nn and pp / 3": "standard rates with nn and pp bremsstrahlung divided by 3 (literature normalization)",
    "updated, no beta": "updated rates without the short-range suppression factors beta",
    "updated, neutron Urca x 1.3": "updated rates with the neutron branch x 1.3 (bound for the non-triangular integral)",
}


def _scaled(function, factor):
    return lambda state: factor * function(state)


def _without_muons(function):
    def wrapped(state):
        stripped = copy.copy(state)
        stripped.kf_mu = np.zeros_like(state.kf_mu)
        return function(stripped)

    return wrapped


def configure(variant):
    """Set the neutrino registry for a variant (starting from the standard rates)."""
    updated.use_standard_rates()
    register = lambda name, function: NEUTRINO_PROCESSES.register(name, replace=True)(
        function
    )  # noqa: E731
    new = updated.make_processes(short_range_suppression=variant != "updated, no beta")
    if variant in ("updated", "updated, no beta"):
        updated.use_updated_rates(short_range_suppression=variant != "updated, no beta")
    elif variant == "updated, no muon channels":
        register("modified_urca", _without_muons(new["modified_urca"]))
        for name in ("brems_nn", "brems_np", "brems_pp"):
            register(name, new[name])
    elif variant == "updated modified Urca only":
        register("modified_urca", new["modified_urca"])
    elif variant == "updated bremsstrahlung only":
        for name in ("brems_nn", "brems_np", "brems_pp"):
            register(name, new[name])
    elif variant == "standard, nn and pp / 3":
        register("brems_nn", _scaled(standard.neutron_neutron_bremsstrahlung, 1.0 / 3.0))
        register("brems_pp", _scaled(standard.proton_proton_bremsstrahlung, 1.0 / 3.0))
    elif variant == "updated, neutron Urca x 1.3":
        updated.use_updated_rates()
        urca = new["modified_urca"]

        def neutron_branch_scaled(state):  # only the neutron branch: total - proton branch
            return urca(state) + 0.3 * _neutron_branch(state)

        register("modified_urca", neutron_branch_scaled)
    elif variant != "standard":
        raise ValueError(f"unknown variant {variant!r}")


def _neutron_branch(state):
    """The neutron branch of the updated modified Urca (with muons), for the sensitivity variant."""
    k = updated._setup(state)  # noqa: SLF001
    f = standard.process_factors(state)
    ve = updated.lepton_velocity(k["p_e"], updated.ELECTRON_MASS_MEV)
    vmu = np.where(k["p_mu"] > 0, updated.lepton_velocity(k["p_mu"], updated.MUON_MASS_MEV), 0.0)
    q = updated.murca_neutron_branch(k["p_n"], k["p_p"], ve + vmu, k["mn"], k["mp"], k["T"])
    out = np.zeros_like(state.T)
    out[k["core"]] = updated.BETA["murca"] * q * f.Rmn[k["core"]] * updated.MEV5_TO_CGS
    return out


def switch_off_superfluidity(text):
    """The text of an input file with the three superfluidity switches (lines a, b, c) set to 0."""
    lines, armed = text.split("\n"), False
    for i, line in enumerate(lines):
        if armed and re.match(r"^\s*1\s+\d+\s*$", line):
            lines[i], armed = "0" + line.lstrip()[1:], False
        if re.match(r"^[abc]\)", line):
            armed = True
    return "\n".join(lines)


@contextmanager
def normal_matter_input(path):
    """A temporary copy of the input file without superfluidity, next to the original (it finds ``eos``)."""
    path = Path(path)
    temporary = path.with_name(f".{path.name}_normal_matter")
    temporary.write_text(switch_off_superfluidity(path.read_text()))
    try:
        yield temporary
    finally:
        temporary.unlink(missing_ok=True)


def luminosity_at(run, ages):
    """L_infinity at the given ages (log-log interpolation; NaN beyond the end of the run)."""
    t, L = np.asarray(run.times_years), np.asarray(run.luminosity_erg_s)
    values = np.exp(np.interp(np.log(ages), np.log(t), np.log(L)))
    return np.where((np.asarray(ages) >= t[0]) & (np.asarray(ages) <= t[-1]), values, np.nan)


def compare(background, final_years=2e5, variants=None):
    """Run every variant; returns ``{variant: CoolingRun}``. The registry is restored afterwards."""
    previous = NEUTRINO_PROCESSES.items()
    runs = {}
    config = replace(CoolingConfig(), final_years=final_years)
    try:
        for variant in variants or VARIANTS:
            configure(variant)
            runs[variant] = run_cooling(background, config)
    finally:
        NEUTRINO_PROCESSES.restore(previous)
    return runs


def write_csv(runs, path):
    with open(path, "w") as f:
        f.write("variant,time_years,L_infinity_erg_s\n")
        for variant, run in runs.items():
            for t, L in zip(run.times_years, run.luminosity_erg_s):
                f.write(f'"{variant}",{t:.10g},{L:.10g}\n')


def plot(runs, path, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    style = {
        "standard": dict(color="black", lw=2.2),
        "updated": dict(color="tab:red", lw=2.2),
        "updated, no muon channels": dict(color="tab:orange", lw=1.4, ls="--"),
        "updated modified Urca only": dict(color="tab:blue", lw=1.3, ls="-."),
        "updated bremsstrahlung only": dict(color="tab:green", lw=1.3, ls="-."),
        "standard, nn and pp / 3": dict(color="gray", lw=1.4, ls=":"),
        "updated, no beta": dict(color="tab:purple", lw=1.3, ls="--"),
        "updated, neutron Urca x 1.3": dict(color="tab:brown", lw=1.3, ls="--"),
    }
    ages = np.logspace(0, 5.3, 400)
    base = luminosity_at(runs["standard"], ages)
    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(8, 8.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 2])
    )
    for variant, run in runs.items():
        top.loglog(run.times_years, run.luminosity_erg_s, label=variant, **style.get(variant, {}))
        bottom.semilogx(
            ages, 100.0 * (luminosity_at(run, ages) / base - 1.0), **style.get(variant, {})
        )
    top.set_ylabel(r"$L_\infty$ [erg s$^{-1}$]")
    top.set_ylim(1e31, 3e35)
    top.legend(fontsize=8, loc="lower left")
    top.set_title(title)
    bottom.axhline(0.0, color="black", lw=0.6)
    bottom.set_xlabel("age [yr]")
    bottom.set_ylabel(r"$L_\infty$ relative to standard [%] (symlog)")
    bottom.set_yscale("symlog", linthresh=10.0)
    peak = max(np.nanmax(100.0 * (luminosity_at(run, ages) / base - 1.0)) for run in runs.values())
    bottom.set_ylim(-30, max(30.0, 1.2 * peak))
    bottom.set_xlim(1, 2e5)
    for ax in (top, bottom):
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def summary(runs, ages=(1, 10, 30, 100, 300, 1e3, 1e4, 1e5)):
    base = luminosity_at(runs["standard"], np.array(ages, float))
    lines = [
        "age [yr]".ljust(32) + "".join(f"{a:>9g}" for a in ages),
        "L_inf standard [erg/s]".ljust(32) + "".join(f"{x:9.2e}" for x in base),
    ]
    for variant, run in runs.items():
        if variant == "standard":
            continue
        ratio = luminosity_at(run, np.array(ages, float)) / base - 1.0
        lines.append((variant + " [%]").ljust(32) + "".join(f"{100 * x:+9.1f}" for x in ratio))
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input", type=Path, default=ROOT / "input_bsk26")
    parser.add_argument(
        "--star", default="1.40 Msun", help="label of the star in the input (default 1.40 Msun)"
    )
    parser.add_argument("--final-years", type=float, default=2e5)
    parser.add_argument(
        "--normal-matter", action="store_true", help="switch the superfluidity of the input off"
    )
    parser.add_argument("--out", type=Path, default=ROOT / "results_emissivities")
    args = parser.parse_args(argv)
    if args.normal_matter:
        with normal_matter_input(args.input) as temporary:
            stars = dict(load_backgrounds(temporary))
    else:
        stars = dict(load_backgrounds(args.input))
    if args.star not in stars:
        parser.error(f"star {args.star!r} not in the input; available: {', '.join(stars)}")
    runs = compare(stars[args.star], args.final_years)
    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(runs, args.out / "emissivity_comparison.csv")
    plot(
        runs,
        args.out / "emissivity_comparison.png",
        f"{args.star}, {args.input.name}{', normal matter' if args.normal_matter else ''}: "
        "standard and updated neutrino emissivities",
    )
    print(summary(runs))
    print(f"Files in {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
