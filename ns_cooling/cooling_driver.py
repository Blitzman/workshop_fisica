"""Run the B=0 cooling of one or several stars, save CSV data and show a synchronized animation."""

import argparse
import csv
import re
from dataclasses import replace
from pathlib import Path

import numpy as np

from background import Background
from config import CoolingConfig
from cooling import run_cooling
from emissivity_updated import use_standard_rates, use_updated_rates
from emissivity import NEUTRINO_PROCESSES, modified_urca_emissivity
from energy_balance import energy_budget, energy_closure, validate_energy_options
from input_model import load_backgrounds, read_parameters

ROOT = Path(__file__).resolve().parent


def folder_name(label):
    """A folder name from the label of a star ("1.40 Msun" -> "1.40_Msun")."""
    return re.sub(r"[^A-Za-z0-9.]+", "_", label).strip("_") or "star"


def save_run(out, background, run):
    """Write the CSV files of one run into the folder ``out``."""
    with open(out / "neutrino_processes.csv", "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_years", "dominant_process"] + ["Lnu_inf_" + name + "_erg_s" for name, _ in run.neutrino_processes])
        for time, dominant, values in zip(run.step_times_years, run.dominant_neutrino_process, run.neutrino_luminosities):
            writer.writerow([time, dominant, *values])
    grid = run.step_times_years
    steps = sorted(run.profiles)
    np.savetxt(
        out / "cooling_curve.csv",
        np.column_stack([run.times_years, run.luminosity_erg_s]),
        delimiter=",",
        header="time_years,L_infinity_erg_s",
        comments="",
    )
    np.savetxt(
        out / "temperature_profiles.csv",
        np.column_stack([background.cells.rho] + [run.profiles[n] * 1e8 for n in steps]),
        delimiter=",",
        header="rho_cell_g_cm3," + ",".join(f"Tinf_K_step{n}" for n in steps),
        comments="",
    )
    np.savetxt(
        out / "profile_times.csv",
        np.column_stack([steps, grid[steps]]),
        delimiter=",",
        header="step,time_years",
        comments="",
        fmt=["%d", "%.16e"],
    )
    np.savetxt(
        out / "time_steps.csv",
        np.column_stack(
            [np.arange(1, run.completed_steps + 1), grid[1:], np.diff(grid), run.relative_changes]
        ),
        delimiter=",",
        header="step,time_years,dt_years,max_relative_temperature_change",
        comments="",
    )


def main():
    defaults = CoolingConfig()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "input_bsk26",
        help="Input file (default: input_bsk26, one 1.4 Msun BSk26 star)",
    )
    parser.add_argument(
        "--min-temperature",
        type=float,
        default=defaults.minimum_temperature_K,
        help="Minimum local interior temperature [K]; the run stops when reached",
    )
    parser.add_argument("--min-luminosity", type=float, default=defaults.minimum_luminosity_erg_s)
    parser.add_argument("--final-years", type=float, default=defaults.final_years)
    parser.add_argument(
        "--predictor-corrector",
        action=argparse.BooleanOptionalAction,
        default=defaults.predictor_corrector,
        help="Predictor-corrector with midpoint coefficients (default on; see README)",
    )
    parser.add_argument(
        "--theta",
        type=float,
        default=defaults.theta,
        help="Theta of the corrector: 0.5 = Crank-Nicolson (default), 1 = backward Euler (see README)",
    )
    parser.add_argument(
        "--background",
        type=Path,
        help="Face table (CSV) to use as the background instead of the structure and "
        "superfluidity lines of the input file (see README, BACKGROUND)",
    )
    parser.add_argument(
        "--save-background",
        type=Path,
        help="Write the background used as a face table (CSV) and continue",
    )
    parser.add_argument(
        "--initial-temperature",
        type=float,
        help="Initial T_infinity in units of 1e8 K; default: the value of the input file",
    )
    parser.add_argument("--out", type=Path, default=ROOT / "results_cooling")
    parser.add_argument(
        "--interval-ms",
        type=float,
        default=40,
        help="Animation interval per frame [ms] (240 frames equally spaced in log t)",
    )
    parser.add_argument("--no-animation", action="store_true", help="Save CSV and final shaded plot without animation")
    parser.add_argument(
        "--check-energy",
        action="store_true",
        help="Check the energy balance dE/dt = -(L_nu + L_gamma) of the run, print the closure "
        "error per decade and write energy_budget.csv",
    )
    parser.add_argument(
        "--energy-start-years",
        type=float,
        default=1.0,
        help="First time of the energy check [yr] (the unresolved start never closes)",
    )
    parser.add_argument(
        "--energy-tolerance",
        type=float,
        default=0.02,
        help="Largest acceptable closure error per decade (fraction of the energy lost)",
    )
    parser.add_argument("--murca", choices=("friman-maxwell", "bottaro"), default="friman-maxwell",
                        help="Modified Urca prescription (default: friman-maxwell); bremsstrahlung unchanged")
    parser.add_argument(
        "--old-rates",
        action="store_true",
        help="Use the standard (Friman-Maxwell) modified Urca and nucleon bremsstrahlung rates "
        "(legacy option selecting all standard rates; incompatible with --murca bottaro)",
    )
    args = parser.parse_args()
    if args.old_rates and args.murca == "bottaro":
        parser.error("--old-rates cannot be combined with --murca bottaro")
    if args.old_rates:
        use_standard_rates()
    else:
        use_updated_rates()
        if args.murca == "friman-maxwell":
            NEUTRINO_PROCESSES.register("modified_urca", replace=True)(modified_urca_emissivity)
    if args.check_energy:
        try:
            validate_energy_options(args.energy_start_years, args.energy_tolerance)
        except ValueError as error:
            parser.error(str(error))
    if not np.isfinite(args.interval_ms) or args.interval_ms <= 0:
        parser.error("--interval-ms must be positive and finite")
    try:
        initial = args.initial_temperature
        if initial is None:
            initial = read_parameters(args.input)["Tinit8"]
        config = replace(
            defaults,
            initial_temperature_8=initial,
            final_years=args.final_years,
            minimum_temperature_K=args.min_temperature,
            minimum_luminosity_erg_s=args.min_luminosity,
            theta=args.theta,
            predictor_corrector=args.predictor_corrector,
        )
    except ValueError as error:
        parser.error(str(error))

    if args.background:
        stars = [("background", Background.from_csv(args.background))]
    else:
        stars = load_backgrounds(args.input)
    several = len(stars) > 1
    out = args.out.resolve()
    entries, curves = [], []
    rates = ("MURCA: Friman-Maxwell; nucleon bremsstrahlung: standard" if args.old_rates
             else f"MURCA: {args.murca}; nucleon bremsstrahlung: Bottaro-Caputo-Fiorillo")
    print(f"Neutrino rates: {rates}")
    for label, background in stars:
        if several:
            print(f"=== {label} ===")
        if args.save_background:
            path = args.save_background
            if several:
                path = path.with_name(f"{path.stem}_{folder_name(label)}{path.suffix}")
            background.to_csv(path)
            print(f"Background written to {path}")
        run = run_cooling(background, config)
        if background.description:
            print(background.description)
        print(f"Envelope: {background.envelope_model}")
        star_out = out / folder_name(label) if several else out
        star_out.mkdir(parents=True, exist_ok=True)
        save_run(star_out, background, run)
        print(f"Accepted steps: {run.completed_steps}; rejected trials: {run.rejected_steps}")
        print(f"Stop: {run.stop_reason}; t={run.step_times_years[run.completed_steps]:.12g} yr")
        print(f"CSV: {star_out}")
        if args.check_energy:
            budget = energy_budget(run, background)
            names, columns = budget.as_columns()
            np.savetxt(
                star_out / "energy_budget.csv",
                columns,
                delimiter=",",
                header=",".join(names),
                comments="",
            )
            print(
                energy_closure(
                    budget, start_years=args.energy_start_years, tolerance=args.energy_tolerance
                ).format()
            )
        entries.append((label, run, background.cells.rho))
        curves.append((label, run))
    if several:
        with open(out / "cooling_curves.csv", "w") as f:
            f.write("star,time_years,L_infinity_erg_s\n")
            for label, run in curves:
                for t, lum in zip(run.times_years, run.luminosity_erg_s):
                    f.write(f"{label},{t:.17g},{lum:.17g}\n")
        print(f"Combined luminosity curves: {out / 'cooling_curves.csv'}")
    from neutrino_dominance import save_dominance_plot
    final_figure = save_dominance_plot(curves, out / "cooling_neutrino_dominance.png", config.minimum_luminosity_erg_s)
    print(f"Final plot: {out / 'cooling_neutrino_dominance.png'}")
    if args.no_animation:
        import matplotlib.pyplot as plt
        plt.close(final_figure)
    if not args.no_animation:
        import matplotlib.pyplot as plt

        from cooling_animation import cooling_animation

        fig, animation = cooling_animation(entries, args.interval_ms)
        plt.show(block=True)


if __name__ == "__main__":
    main()
