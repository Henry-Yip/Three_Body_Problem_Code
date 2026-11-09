"""
Makes the figures in the README and saves them in the figures folder.

    python make_figures.py

Author: Henry Yip
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from fixed_centres import ORBITS, conserved, simulate
from three_body import equal_mass_start, integrate, total_energy

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
COLOURS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")

# One colour per orbit type, used in every figure.
SATELLITE, LEMNISCATE, PLANETARY = COLOURS[:3]

# p1, p2 and the period of each orbit, from Suvakov and Dmitrasinovic.
PERIODIC_ORBITS = {
    "Figure-8": (0.347111, 0.532728, 6.3259),
    "Moth I": (0.464445, 0.396060, 14.8939),
    "Bumblebee": (0.184279, 0.587188, 63.5345),
}

plt.rcParams.update({
    "font.family": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.labelcolor": SECONDARY,
    "axes.titlecolor": INK,
    "axes.titlesize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "text.color": INK,
    "xtick.color": AXIS,
    "ytick.color": AXIS,
    "xtick.labelcolor": SECONDARY,
    "ytick.labelcolor": SECONDARY,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "savefig.dpi": 200,
    "savefig.facecolor": SURFACE,
    "lines.solid_capstyle": "round",
    "lines.solid_joinstyle": "round",
})


def periodic_orbits(folder):
    """The paths of the three bodies over one period of three of the orbits in Table 1."""
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.4))
    dt = 1e-4
    every = 10
    for ax, (name, (p1, p2, period)) in zip(axes, PERIODIC_ORBITS.items()):
        steps = round(period / dt)
        positions, velocities = equal_mass_start(p1, p2)
        path = np.zeros((steps // every + 1, 3, 2))
        for i, (x, _) in enumerate(integrate(positions, velocities, np.ones(3), dt, steps)):
            if i % every == 0:
                path[i // every] = x[:, :2]
        path = path[: steps // every]
        for body in range(3):
            ax.plot(path[:, body, 0], path[:, body, 1], color=COLOURS[body], linewidth=0.9)
            ax.plot(*np.array(positions)[body, :2], "o", color=COLOURS[body], markersize=5,
                    markeredgecolor=SURFACE, markeredgewidth=1)
        ax.set_title(f"{name}, period {period:g}")
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    keys = [Line2D([], [], color=COLOURS[body], marker="o", markersize=5, linewidth=1.5,
                   markeredgecolor=SURFACE, label=f"Body {body + 1}") for body in range(3)]
    fig.legend(handles=keys, loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(folder, "periodic_orbits.png"))
    plt.close(fig)


def convergence(folder):
    """How the energy error of each method falls with the step size, on the Moth I orbit."""
    methods = (("neri", "Neri", 4), ("ruth", "Ruth", 3), ("verlet", "Verlet", 2), ("euler", "Euler", 1))
    step_sizes = np.array([0.004, 0.002, 0.001, 0.0005, 0.00025, 0.000125])
    p1, p2, _ = PERIODIC_ORBITS["Moth I"]
    masses = np.ones(3)
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for colour, (method, name, order) in zip(COLOURS, methods):
        errors = []
        for dt in step_sizes:
            positions, velocities = equal_mass_start(p1, p2)
            orbit = integrate(positions, velocities, masses, dt, round(10 / dt), method)
            energy = np.array([total_energy(x, v, masses) for x, v in orbit])
            errors.append(abs((energy.max() - energy.min()) / energy[0]))
        ax.plot(step_sizes, errors, color=colour, linewidth=2, marker="o", markersize=6,
                markeredgecolor=SURFACE, markeredgewidth=1.2, label=name)
        ax.annotate(f"order {order}", (step_sizes[-1], errors[-1]), xytext=(-8, 0),
                    textcoords="offset points", va="center", ha="right", color=SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(step_sizes)
    ax.set_xticklabels([f"{dt:g}" for dt in step_sizes])
    ax.minorticks_off()
    ax.set_xlim(step_sizes[-1] * 0.5, step_sizes[0] * 1.3)
    ax.set_xlabel("Step size")
    ax.set_ylabel("Energy deviation")
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", ncol=1)
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "convergence.png"))
    plt.close(fig)


def orbit_types(folder, values):
    """The regions of the (h, g) plane that hold each orbit type, for unit masses and b = 1."""
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    h = np.linspace(-2, -1e-3, 2000)
    top = 4.0
    planetary = h >= -1
    ax.fill_between(h, h, 0, color=SATELLITE, alpha=0.14, linewidth=0)
    ax.fill_between(h, 0, 2 + h, color=LEMNISCATE, alpha=0.14, linewidth=0)
    ax.fill_between(h[planetary], 2 + h[planetary], np.minimum(1 / np.abs(h[planetary]), top),
                    color=PLANETARY, alpha=0.14, linewidth=0)

    # The solid lines are the boundaries that follow from the separated equations.
    ax.plot(h, h, color=SECONDARY, linewidth=1)
    ax.plot(h, 0 * h, color=SECONDARY, linewidth=1)
    ax.plot(h, 2 + h, color=SECONDARY, linewidth=1)
    outer = 1 / np.abs(h[planetary])
    ax.plot(h[planetary][outer <= top], outer[outer <= top], color=SECONDARY, linewidth=1)
    ax.plot([-1, 0], [1, 1], color=INK, linewidth=1.2, linestyle=(0, (4, 3)))

    ax.text(-1.25, -0.5, "Satellite", ha="center", va="center")
    ax.text(-1.05, 0.42, "Lemniscate", ha="center", va="center")
    ax.text(-0.22, 2.75, "Planetary", ha="center", va="center")
    ax.text(-1.45, 2.4, "No orbits", ha="center", va="center", color=MUTED)
    ax.text(-1.6, 0.5, r"$g = 2 - |h|$", color=SECONDARY, ha="center", va="bottom", rotation=45,
            rotation_mode="anchor", transform_rotates_text=True)
    ax.text(-0.03, 1.06, "the paper's cut at g = 1", color=SECONDARY, ha="right", va="bottom")

    colours = {"planetary": PLANETARY, "lemniscate": LEMNISCATE, "satellite": SATELLITE}
    for name, (g, h_value) in values.items():
        ax.plot(h_value, g, "o", color=colours[name], markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.5)

    ax.set_xlim(-2, 0)
    ax.set_ylim(-2, top)
    ax.set_xlabel("h (energy)")
    ax.set_ylabel("g (Euler's second integral)")
    ax.grid(True)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "orbit_types.png"))
    plt.close(fig)


def planetary_orbit(folder, x, z, dt, g, h):
    """The planetary orbit in the plane, and R over time."""
    a = 1 / abs(h)
    e = np.sqrt(1 - g / a)
    r_min, r_max = a * (1 - e), a * (1 + e)
    fig, (left, right) = plt.subplots(1, 2, figsize=(9, 4.2), gridspec_kw={"width_ratios": [1, 1.25]})

    left.plot(x, z, color=PLANETARY, linewidth=0.8)
    angle = np.linspace(0, 2 * np.pi, 400)
    for R in (r_min, r_max):
        left.plot(np.sqrt(R * R - 1) * np.sin(angle), R * np.cos(angle), color=SECONDARY, linewidth=1.2)
    left.plot([0, 0], [1, -1], "o", color=INK, markersize=5, markeredgecolor=SURFACE, markeredgewidth=1)
    left.text(0, 0, "Fixed\ncentres", ha="center", va="center", color=SECONDARY)
    left.text(0, r_min - 0.12, f"R = {r_min:.2f}", ha="center", va="top", color=SECONDARY)
    left.text(0, r_max + 0.1, f"R = {r_max:.2f}", ha="center", va="bottom", color=SECONDARY)
    left.set_ylim(-3.7, 4.0)
    left.set_aspect("equal")
    left.set_xlabel("x")
    left.set_ylabel("z")
    left.set_title("Planetary orbit in the x–z plane")

    time = dt * (np.arange(len(x)) + 1)
    R = 0.5 * (np.hypot(x, z - 1) + np.hypot(x, z + 1))
    right.plot(time, R, color=PLANETARY, linewidth=1.5)
    right.axhline(a, color=INK, linewidth=1)
    right.axhline(g, color=INK, linewidth=1, linestyle=(0, (4, 3)))
    end = time[-1]
    right.text(end * 1.015, a, rf"$a = 1/|h| = {a:.3f}$", va="center", color=SECONDARY)
    right.text(end * 1.015, g, rf"$p = {g:.3f}$", va="center", color=SECONDARY)
    right.set_xlim(0, end)
    right.set_xlabel("Time")
    right.set_ylabel("R")
    right.grid(True, axis="y")
    right.set_axisbelow(True)
    right.set_title("R against time")
    fig.tight_layout()
    fig.savefig(os.path.join(folder, "planetary_orbit.png"))
    plt.close(fig)


def main():
    folder = "figures"
    os.makedirs(folder, exist_ok=True)

    values = {}
    for name, (start, velocity, dt, steps) in ORBITS.items():
        x, z, vx, vz = simulate(start, velocity, dt, steps)
        g, h = conserved(x, z, vx, vz)
        values[name] = (g.mean(), h.mean())
        if name == "planetary":
            planetary_orbit(folder, x, z, dt, *values[name])

    orbit_types(folder, values)
    convergence(folder)
    periodic_orbits(folder)


if __name__ == "__main__":
    main()
