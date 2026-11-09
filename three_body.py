"""
Symplectic integration of point masses under Newtonian gravity, with G = 1.

The scheme is the one set out in H. Yip and J. M. Smillie, "An Analysis of the
Restricted Euler Problem Using Symplectic Integrators", Edinburgh Student Journal
of Science 1(3), 2025. Each step runs through a list of stages, and each stage does
    v = v + c * a(x) * dt
    x = x + d * v * dt

To save a trajectory as an XYZ file:
    python three_body.py 0.347111 0.532728 100000 0.0001 figure_8.xyz --every 10
The two numbers at the start are the velocity (p1, p2) of the outer bodies, which
start at (1, 0) and (-1, 0). The third body starts at the origin with velocity
(-2 p1, -2 p2), so the total momentum is zero.

Author: Henry Yip
"""
import argparse

import numpy as np

_CBRT2 = 2 ** (1 / 3)
_W = 2 - _CBRT2

# (c coefficients, d coefficients) for each method.
METHODS = {
    "neri": ([1 / (2 * _W), (1 - _CBRT2) / (2 * _W), (1 - _CBRT2) / (2 * _W), 1 / (2 * _W)],
             [1 / _W, -_CBRT2 / _W, 1 / _W, 0.0]),
    "ruth": ([7 / 24, 3 / 4, -1 / 24], [2 / 3, -2 / 3, 1.0]),
    "verlet": ([1 / 2, 1 / 2], [1.0, 0.0]),
    # Position first, then velocity. This is the order that produced Table 1.
    "euler": ([0.0, 1.0], [1.0, 0.0]),
}


def accelerations(positions, masses):
    """Return the acceleration of every body due to all the other bodies."""
    separation = positions[np.newaxis, :, :] - positions[:, np.newaxis, :]
    distance_squared = np.sum(separation ** 2, axis=2)
    np.fill_diagonal(distance_squared, np.inf)
    weight = masses[np.newaxis, :] / distance_squared ** 1.5
    return np.sum(weight[:, :, np.newaxis] * separation, axis=1)


def total_energy(positions, velocities, masses):
    """Return the kinetic energy plus the gravitational potential energy."""
    kinetic = 0.5 * np.sum(masses * np.sum(velocities ** 2, axis=1))
    first, second = np.triu_indices(len(masses), k=1)
    distances = np.linalg.norm(positions[first] - positions[second], axis=1)
    return kinetic - np.sum(masses[first] * masses[second] / distances)


def integrate(positions, velocities, masses, dt, steps, method="neri", free=None):
    """
    Yield (positions, velocities) after each step.

    The optional array `free` marks the bodies that are allowed to move. Any other
    body stays where it is and still attracts the others, which is how the two fixed
    centres are modelled. The arrays yielded are updated in place between steps.
    """
    c_list, d_list = METHODS[method]
    positions = np.array(positions, dtype=float)
    velocities = np.array(velocities, dtype=float)
    masses = np.asarray(masses, dtype=float)
    free = np.ones(len(masses), dtype=bool) if free is None else np.asarray(free, dtype=bool)

    for _ in range(steps):
        for c, d in zip(c_list, d_list):
            if c != 0:
                velocities[free] += c * dt * accelerations(positions, masses)[free]
            if d != 0:
                positions[free] += d * dt * velocities[free]
        yield positions, velocities


def write_frame(handle, number, positions, labels):
    """Write one frame in the XYZ format described in the README."""
    handle.write(f"{len(labels)}\nPoint = {number}\n")
    for label, (x, y, z) in zip(labels, positions):
        handle.write(f"{label}    {x} {y} {z}\n")


def equal_mass_start(p1, p2):
    """Return the starting positions and velocities used for the periodic orbits."""
    positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
    velocities = [[-2 * p1, -2 * p2, 0.0], [p1, p2, 0.0], [p1, p2, 0.0]]
    return positions, velocities


def main():
    parser = argparse.ArgumentParser(description="Save an equal-mass three-body orbit as an XYZ file.")
    parser.add_argument("p1", type=float)
    parser.add_argument("p2", type=float)
    parser.add_argument("steps", type=int)
    parser.add_argument("dt", type=float)
    parser.add_argument("output")
    parser.add_argument("--method", choices=METHODS, default="neri")
    parser.add_argument("--every", type=int, default=1, help="save every Nth step")
    args = parser.parse_args()

    positions, velocities = equal_mass_start(args.p1, args.p2)
    labels = ("Central_Planet", "First_Planet", "Second_Planet")
    orbit = integrate(positions, velocities, np.ones(3), args.dt, args.steps, args.method)
    with open(args.output, "w") as handle:
        for step, (x, _) in enumerate(orbit):
            if step % args.every == 0:
                write_frame(handle, step // args.every + 1, x, labels)


if __name__ == "__main__":
    main()
