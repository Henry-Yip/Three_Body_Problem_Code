"""
The problem of two fixed centres: the three orbits used in the paper, their conserved
quantities, and the measurement of the planetary orbit's semi-major axis.

    python fixed_centres.py          prints the results
    python fixed_centres.py --xyz    also saves each trajectory as an XYZ file

The centres have unit mass and sit at (x, z) = (0, 1) and (0, -1), with G = 1.

Author: Henry Yip
"""
import argparse

import numpy as np

from three_body import integrate, write_frame

MASSES = np.ones(3)
MOVING = np.array([True, False, False])
CENTRES = [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]
LABELS = ("Central_Planet", "First_Planet", "Second_Planet")
SAVE_INTERVAL = 0.01

# Starting position and velocity (x, y, z) of the moving body, the step size and the
# number of steps. The lemniscate and satellite orbits pass close to a centre, so they
# need the smaller step.
ORBITS = {
    "planetary": ((1.8, 0.0, 1.0), (2 / 3, 0.0, -0.8), 1e-3, 120_000),
    "lemniscate": ((2.5, 0.0, 1.5), (0.0, 0.0, -0.4), 1e-4, 400_000),
    "satellite": ((2.5, 0.0, 2.0), (0.0, 0.0, 0.0), 1e-4, 400_000),
}


def conserved(x, z, vx, vz):
    """Return the paper's g and h for b = 1 and unit masses."""
    to_lower = np.hypot(x, z + 1)
    to_upper = np.hypot(x, z - 1)
    g = 0.5 * (z * vx - x * vz) ** 2 + 0.5 * vz ** 2 + z * (1 / to_lower - 1 / to_upper)
    h = 0.5 * (vx ** 2 + vz ** 2) - 1 / to_lower - 1 / to_upper
    return g, h


def simulate(start, velocity, dt, steps):
    """Return x, z, vx and vz of the moving body after each step."""
    positions = [start, *CENTRES]
    velocities = [velocity, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    path = np.zeros((steps, 4))
    for i, (r, v) in enumerate(integrate(positions, velocities, MASSES, dt, steps, "neri", MOVING)):
        path[i] = r[0, 0], r[0, 2], v[0, 0], v[0, 2]
    return path.T


def orbit_type(x, z):
    """Classify a trajectory by where it crosses the line through the two centres."""
    if np.all(z > 0) or np.all(z < 0):
        return "satellite"
    i = np.nonzero(x[:-1] * x[1:] < 0)[0]
    crossing_z = z[i] + (z[i + 1] - z[i]) * x[i] / (x[i] - x[i + 1])
    return "lemniscate" if np.any(np.abs(crossing_z) < 1) else "planetary"


def half_extents(z, loops=6):
    """Return half the distance between the top and bottom of each of the first loops."""
    inner = z[1:-1]
    peaks = np.nonzero((inner > z[:-2]) & (inner >= z[2:]) & (inner > 0))[0] + 1
    troughs = np.nonzero((inner < z[:-2]) & (inner <= z[2:]) & (inner < 0))[0] + 1
    count = min(len(peaks), len(troughs), loops)
    return 0.5 * (z[peaks[:count]] - z[troughs[:count]])


def save_xyz(path, x, z, dt):
    every = round(SAVE_INTERVAL / dt)
    with open(path, "w") as handle:
        for frame, step in enumerate(range(0, len(x), every), start=1):
            write_frame(handle, frame, [[x[step], 0.0, z[step]], *CENTRES], LABELS)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--xyz", action="store_true", help="save each trajectory as an XYZ file")
    args = parser.parse_args()

    print(f"{'Orbit':12}{'g':>9}{'h':>9}{'spread g':>11}{'spread h':>11}{'p = |g|':>9}{'eta':>8}{'2 - |h|':>9}  Type from trajectory")
    for name, (start, velocity, dt, steps) in ORBITS.items():
        x, z, vx, vz = simulate(start, velocity, dt, steps)
        g_values, h_values = conserved(x, z, vx, vz)
        g, h = g_values.mean(), h_values.mean()
        print(f"{name:12}{g:9.4f}{h:9.4f}{np.ptp(g_values):11.1e}{np.ptp(h_values):11.1e}"
              f"{abs(g):9.4f}{1 / abs(g):8.4f}{2 - abs(h):9.4f}  {orbit_type(x, z)}")

        if name == "planetary":
            semi_major = 1 / abs(h)
            eccentricity = np.sqrt(1 - g / semi_major)
            mean_distance = 0.5 * (np.hypot(x, z - 1) + np.hypot(x, z + 1))
            loops = half_extents(z)
            planetary = (semi_major, eccentricity, mean_distance.min(), mean_distance.max(), loops)

        if args.xyz:
            save_xyz(f"{name}.xyz", x, z, dt)
            print(f"    saved {name}.xyz, with {SAVE_INTERVAL:g} between saved points")

    a, e, r_min, r_max, loops = planetary
    print("\nPlanetary orbit")
    print(f"  semi-major axis a = 1/|h| = {a:.4f}, and eccentricity e = sqrt(1 - p/a) = {e:.4f}")
    print(f"  R runs from {r_min:.3f} to {r_max:.3f}; the book's values a(1 - e) and a(1 + e) are {a * (1 - e):.3f} and {a * (1 + e):.3f}")
    print(f"  half-extent in z of the first six loops: {np.round(loops, 2)}")
    print(f"  the mean of the six is {loops.mean():.2f}, and the population standard deviation of the last five is {loops[1:].std():.3f}")
    print("  the paper reports a mean of 2.53 with a standard deviation of 0.133")


if __name__ == "__main__":
    main()
