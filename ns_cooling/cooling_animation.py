"""Synchronized playback of computed cooling histories (no solver changes).

Several stars can be compared. Each run has its own adaptive time steps, so the
playback uses a common grid of times, equally spaced in log t, and shows every star
at the same physical time: the luminosity curve up to that time and the redshifted
temperature profile at that time (interpolated linearly in log t and log T between the
accepted steps). A star whose run ended earlier keeps its final state.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button

COLORS = ("tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown")
DEFAULT_FRAMES = 240


def profile_log_temperature(run, times):
    """log10 of T_inf [K] in every cell at the given times (shape (len(times), n_cells)).

    Linear in log t and log T between the accepted steps (the initial state at t = 0 has no
    place on a log axis and is skipped); times beyond the end of the run give the final state.
    """
    steps = sorted(run.profiles)
    t = run.step_times_years[steps]
    keep = t > 0
    log_t = np.log(t[keep])
    log_T = np.log(np.array([run.profiles[n] for n, k in zip(steps, keep) if k]) * 1e8)
    x = np.log(np.minimum(np.asarray(times, dtype=float), t[-1]))
    j = np.clip(np.searchsorted(log_t, x) - 1, 0, len(log_t) - 2)
    w = np.clip((x - log_t[j]) / (log_t[j + 1] - log_t[j]), 0.0, 1.0)[:, None]
    return ((1.0 - w) * log_T[j] + w * log_T[j + 1]) / np.log(10.0)


def frame_data(entries, frames=DEFAULT_FRAMES):
    """Common frame times and what every star shows in each frame.

    ``entries`` is a list of ``(label, run, density)`` with the density of the cells.
    Returns ``(times, stars)``: ``times`` are ``frames`` values equally spaced in log t from the
    first accepted step of any run to the last time of any run; ``stars`` is a list of dicts
    with ``label``, ``run``, ``log_rho`` (cells sorted by density), ``log_T`` (frames x cells,
    same order) and ``count`` (number of luminosity samples reached at each frame).
    """
    if not entries:
        raise ValueError("There is nothing to animate")
    if frames < 2:
        raise ValueError("At least two frames are needed")
    first = min(run.step_times_years[run.step_times_years > 0].min() for _, run, _ in entries)
    last = max(run.step_times_years[run.completed_steps] for _, run, _ in entries)
    times = np.geomspace(first, max(last, first * 1.0001), frames)
    stars = []
    for label, run, density in entries:
        order = np.argsort(density)
        stars.append(
            dict(
                label=label,
                run=run,
                log_rho=np.log10(np.asarray(density, dtype=float)[order]),
                log_T=profile_log_temperature(run, times)[:, order],
                count=np.searchsorted(run.times_years, times, side="right"),
            )
        )
    return times, stars


def cooling_animation(entries, interval_ms=40, frames=DEFAULT_FRAMES):
    """Animate L at infinity and the redshifted temperature profile of one or several stars.

    ``entries`` is a list of ``(label, run, density)`` (see ``frame_data``). Returns
    ``(figure, animation)``.
    """
    times, stars = frame_data(entries, frames)
    fig, (ax_l, ax_t) = plt.subplots(1, 2, figsize=(12.5, 5.8))
    fig.subplots_adjust(bottom=0.20, top=0.85, wspace=0.30)
    title = fig.suptitle("")
    ax_l.set(
        xscale="log",
        yscale="log",
        xlabel="Time [yr]",
        ylabel=r"$L^\infty$ [erg s$^{-1}$]",
        title="Cooling curve",
    )
    sampled = [s["run"] for s in stars if len(s["run"].times_years)]
    if sampled:
        lo = min(r.times_years[0] for r in sampled)
        hi = max(r.times_years[-1] for r in sampled)
        ax_l.set_xlim(lo / 1.1, max(hi * 1.1, lo * 1.2))
        ax_l.set_ylim(
            min(r.luminosity_erg_s.min() for r in sampled) / 1.3,
            max(r.luminosity_erg_s.max() for r in sampled) * 1.3,
        )
    else:
        ax_l.set_xlim(0.01, 1)
        ax_l.set_ylim(1e30, 1e40)
        ax_l.text(0.5, 0.5, "No luminosity samples", transform=ax_l.transAxes, ha="center")
    all_rho = np.concatenate([s["log_rho"] for s in stars])
    all_T = np.concatenate([s["log_T"].ravel() for s in stars])
    ax_t.set(
        xlabel=r"$\log_{10}(\rho\,/\,\mathrm{g\,cm^{-3}})$",
        ylabel=r"$\log_{10}(T^\infty\,/\,\mathrm{K})$",
        title="Temperature profile",
        xlim=(all_rho.min(), all_rho.max()),
        ylim=(all_T.min() - 0.08, all_T.max() + 0.08),
    )
    for k, s in enumerate(stars):
        color = COLORS[k % len(COLORS)]
        (s["curve"],) = ax_l.plot([], [], color=color, lw=2, label=s["label"])
        (s["point"],) = ax_l.plot([], [], "o", color=color, ms=5)
        (s["profile"],) = ax_t.plot([], [], color=color, lw=2, label=s["label"])
    if len(stars) > 1:
        ax_l.legend(loc="upper right", fontsize=9)
    for ax in (ax_l, ax_t):
        ax.grid(alpha=0.25, which="both")
    state = {"index": 0, "paused": False}

    def draw(i):
        artists = []
        for s in stars:
            run, count = s["run"], int(s["count"][i])
            s["curve"].set_data(run.times_years[:count], run.luminosity_erg_s[:count])
            s["point"].set_data(
                run.times_years[count - 1 : count] if count else [],
                run.luminosity_erg_s[count - 1 : count] if count else [],
            )
            s["profile"].set_data(s["log_rho"], s["log_T"][i])
            artists += [s["curve"], s["point"], s["profile"]]
        status = " | Final" if i == len(times) - 1 else ""
        title.set_text(f"Frame {i + 1}/{len(times)} | t = {times[i]:.6g} yr{status}")
        return (*artists, title)

    pause = Button(fig.add_axes([0.38, 0.045, 0.11, 0.055]), "Pause")
    restart = Button(fig.add_axes([0.51, 0.045, 0.11, 0.055]), "Restart")

    def tick(_):
        if not state["paused"] and state["index"] < len(times) - 1:
            state["index"] += 1
        if state["index"] == len(times) - 1:
            animation.event_source.stop()
        return draw(state["index"])

    animation = FuncAnimation(
        fig,
        tick,
        init_func=lambda: draw(0),
        interval=interval_ms,
        blit=False,
        cache_frame_data=False,
    )

    def toggle(_):
        state["paused"] = not state["paused"]
        pause.label.set_text("Resume" if state["paused"] else "Pause")
        fig.canvas.draw_idle()

    def reset(_):
        state.update(index=0, paused=False)
        pause.label.set_text("Pause")
        draw(0)
        animation.event_source.start()
        fig.canvas.draw_idle()

    pause.on_clicked(toggle)
    restart.on_clicked(reset)
    # Keep animation and widgets alive for the lifetime of the window.
    fig._cooling_controls = animation, pause, restart, draw
    draw(0)
    return fig, animation
