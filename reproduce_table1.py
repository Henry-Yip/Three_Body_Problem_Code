"""
Reproduces the energy deviations and momentum differences in Table 1 of the paper.

Each of the twelve runs has 100,000 steps of size 0.0001. The deviation is
|E_max - E_min| / |E_0|, and as in the original runs E_0 is the energy after the
first step. The "Moth" orbit in the paper is Moth I of Suvakov and Dmitrasinovic.

Author: Henry Yip
"""
import numpy as np

from three_body import equal_mass_start, integrate, total_energy

ORBITS = {
    "Figure-8": (0.347111, 0.532728),
    "Bumblebee": (0.184279, 0.587188),
    "Moth": (0.464445, 0.396060),
}

# Energy deviations as printed in the paper.
PAPER = {
    ("neri", "Figure-8"): 6.124e-14,
    ("ruth", "Figure-8"): 1.280e-13,
    ("verlet", "Figure-8"): 5.893e-9,
    ("euler", "Figure-8"): 3.601e-5,
    ("neri", "Bumblebee"): 8.058e-3,
    ("ruth", "Bumblebee"): 8.058e-3,
    ("verlet", "Bumblebee"): 9.676e-2,
    ("euler", "Bumblebee"): 2.4329,
    ("neri", "Moth"): 1.459e-9,
    ("ruth", "Moth"): 2.143e-8,
    ("verlet", "Moth"): 4.775e-5,
    ("euler", "Moth"): 0.02070,
}

STEPS = 100_000
DT = 1e-4


def run(method, p1, p2):
    """Return the energy deviation and the largest changes in x and y momentum."""
    positions, velocities = equal_mass_start(p1, p2)
    masses = np.ones(3)
    energy = np.zeros(STEPS)
    momentum_x = np.zeros(STEPS)
    momentum_y = np.zeros(STEPS)
    for i, (x, v) in enumerate(integrate(positions, velocities, masses, DT, STEPS, method)):
        energy[i] = total_energy(x, v, masses)
        momentum_x[i] = np.sum(v[:, 0])
        momentum_y[i] = np.sum(v[:, 1])
    deviation = abs((energy.max() - energy.min()) / energy[0])
    return deviation, np.ptp(momentum_x), np.ptp(momentum_y)


def main():
    print(f"{'Method':8}{'Orbit':11}{'Energy dev.':>13}{'Paper':>12}{'Max dpx':>12}{'Max dpy':>12}  Comparison")
    for orbit, (p1, p2) in ORBITS.items():
        for method in ("neri", "ruth", "verlet", "euler"):
            deviation, dpx, dpy = run(method, p1, p2)
            printed = PAPER[(method, orbit)]
            if deviation < 1e-12:
                comparison = "round-off level"
            elif abs(deviation / printed - 1) < 1e-3:
                comparison = "matches"
            else:
                comparison = "DIFFERS, see README"
            print(f"{method:8}{orbit:11}{deviation:13.3e}{printed:12.3e}{dpx:12.3e}{dpy:12.3e}  {comparison}", flush=True)


if __name__ == "__main__":
    main()
