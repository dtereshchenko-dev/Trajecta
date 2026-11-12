"""Compare Holman state and STM with ASSIST using existing DE440/SB441 files."""
import argparse
import json
import math
import os
import subprocess

import assist
import numpy as np
import rebound

parser = argparse.ArgumentParser()
parser.add_argument("trajecta", help="Path to the Trajecta comparison driver")
parser.add_argument("--check-fd", action="store_true", help="Check STM against perturbed trajectories")
parser.add_argument("--ias15", action="store_true", help="Use Trajecta IAS15 instead of RA15")
args = parser.parse_args()
de440, sb441 = os.environ["DE440_FILE"], os.environ["SB441_FILE"]
ephem = assist.Ephem(de440, sb441)
ephem.jd_ref = 2458849.5
sun = ephem.get_particle("Sun", 0.0)
initial = np.array([3.338875348598862, -0.9176518412197102, -0.5038590741719294,
                    0.002805663364339457, 0.007550408665778840, 0.002980028207875623])
obliquity = math.radians(84381.448 / 3600.0)
c, s = math.cos(obliquity), math.sin(obliquity)
ecl_to_eq = np.array([[1., 0., 0.], [0., c, -s], [0., s, c]])
eq_to_ecl = ecl_to_eq.T
rotation = np.zeros((6, 6))
rotation[:3, :3] = eq_to_ecl
rotation[3:, 3:] = eq_to_ecl
sun_state = np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz])
initial_ecl = rotation @ (initial - sun_state)
command = [args.trajecta, de440, sb441] + [format(x, ".17g") for x in initial_ecl] + (["--ias15"] if args.ias15 else [])
out = subprocess.check_output(command, text=True)
rows = np.array([[float(value) for value in line.split()] for line in out.splitlines()])
if rows.shape != (6, 43):
    raise RuntimeError(f"unexpected Trajecta output shape {rows.shape}")

sim = rebound.Simulation()
sim.add(x=initial[0], y=initial[1], z=initial[2], vx=initial[3], vy=initial[4], vz=initial[5])
variations = []
names = ("x", "y", "z", "vx", "vy", "vz")
for column in range(6):
    var = sim.add_variation(testparticle=0)
    for row, name in enumerate(names):
        setattr(var.particles[0], name, rotation.T[row, column])
    variations.append(var)
extras = assist.Extras(sim, ephem)
extras.gr_eih_sources = 11
extras.forces = [f for f in extras.forces if f not in ("EARTH_HARMONICS", "NON_GRAVITATIONAL")]

output = {"assist_version": assist.__version__, "rebound_version": rebound.__version__,
          "forces": extras.forces, "gr_eih_sources": extras.gr_eih_sources,
          "epsilon": sim.ri_ias15.epsilon, "min_dt_days": sim.ri_ias15.min_dt,
          "adaptive_mode": str(sim.ri_ias15.adaptive_mode), "results": []}
by_day = {row[0]: row for row in rows}
sample_days = sorted(set(by_day) | set(np.linspace(0., 10000., 501)))
minimum_earth_distance = float("inf")
for day in sample_days:
    extras.integrate_or_interpolate(day)
    p = sim.particles[0]
    earth = ephem.get_particle("Earth", day)
    minimum_earth_distance = min(minimum_earth_distance,
                                 np.linalg.norm(np.array([p.x-earth.x, p.y-earth.y, p.z-earth.z])))
    if day not in by_day:
        continue
    row = by_day[day]
    sun = ephem.get_particle("Sun", day)
    state = rotation @ (np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz]) -
                        np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz]))
    stm = np.column_stack([rotation @ np.array([getattr(v.particles[0], n) for n in names])
                           for v in variations])
    delta = row[1:7] - state
    delta_stm = row[7:].reshape(6, 6) - stm
    output["results"].append({"day": day,
        "position_difference_m": float(np.linalg.norm(delta[:3]) * 149597870700),
        "velocity_difference_m_s": float(np.linalg.norm(delta[3:]) * 149597870700 / 86400),
        "stm_max_abs": float(np.max(np.abs(delta_stm))),
        "stm_position_response_m_for_1e-8_au": float(np.max(np.linalg.norm(delta_stm[:3, :3], axis=0)) * 1e-8 * 149597870700),
        "stm_position_response_m_for_1e-10_au_day": float(np.max(np.linalg.norm(delta_stm[:3, 3:], axis=0)) * 1e-10 * 149597870700),
        "stm_velocity_response_m_s_for_1e-8_au": float(np.max(np.linalg.norm(delta_stm[3:, :3], axis=0)) * 1e-8 * 149597870700 / 86400),
        "stm_velocity_response_m_s_for_1e-10_au_day": float(np.max(np.linalg.norm(delta_stm[3:, 3:], axis=0)) * 1e-10 * 149597870700 / 86400)})
output["minimum_sampled_earth_distance_au"] = float(minimum_earth_distance)
output["earth_distance_sample_count"] = len(sample_days)

if args.check_fd:
    def assist_final(perturbation):
        perturbed = initial + rotation.T @ perturbation
        check = rebound.Simulation()
        check.add(x=perturbed[0], y=perturbed[1], z=perturbed[2],
                  vx=perturbed[3], vy=perturbed[4], vz=perturbed[5])
        check_extra = assist.Extras(check, ephem)
        check_extra.gr_eih_sources = 11
        check_extra.forces = extras.forces
        check_extra.integrate_or_interpolate(10000.)
        p = check.particles[0]
        sun = ephem.get_particle("Sun", 10000.)
        return rotation @ (np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz]) -
                           np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz]))

    def trajecta_final(perturbation):
        cmd = [args.trajecta, de440, sb441] + [format(x, ".17g") for x in initial_ecl + perturbation] + (["--ias15"] if args.ias15 else [])
        raw = subprocess.check_output(cmd, text=True)
        return np.array([float(v) for v in raw.splitlines()[-1].split()])[1:7]

    trajectory_stm = rows[-1, 7:].reshape(6, 6)
    output["finite_difference_checks"] = []
    for column, h in ((0, 1e-7), (0, 1e-8), (3, 1e-9), (3, 1e-10)):
        perturbation = np.zeros(6)
        perturbation[column] = h
        check = {"column": names[column], "perturbation": h}
        for name, propagate, matrix in (("trajecta", trajecta_final, trajectory_stm),
                                         ("assist", assist_final, stm)):
            displacement = (propagate(perturbation) - propagate(-perturbation)) / 2
            check[name + "_position_residual_m"] = float(
                np.linalg.norm(displacement[:3] - matrix[:3, column] * h) * 149597870700)
        output["finite_difference_checks"].append(check)
print(json.dumps(output, indent=2))
