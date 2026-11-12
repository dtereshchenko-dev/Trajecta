"""Compare Parijskij's Ceres encounter with the ASSIST example's sample grid."""
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
parser.add_argument("--convergence", action="store_true", help="Run fixed-step and IAS15 epsilon checks")
parser.add_argument("--output-schedule-check", action="store_true", help="Compare final-only and sampled output")
parser.add_argument("--ias15-check", action="store_true", help="Compare Trajecta IAS15 modes and ASSIST")
parser.add_argument("--ra15-check", action="store_true", help="Compare adaptive RA15 initial steps")
args = parser.parse_args()
de440, sb441 = os.environ["DE440_FILE"], os.environ["SB441_FILE"]
ephem = assist.Ephem(de440, sb441)
ephem.jd_ref = 2449718.5
sun = ephem.get_particle("Sun", 0.0)
initial = np.array([-2.232847879711731, 1.574146331186095, 0.8329414259670296,
                    -0.006247432571575564, -0.007431073424167182, -0.003231725223736132])
obliquity = math.radians(84381.448 / 3600.0)
c, s = math.cos(obliquity), math.sin(obliquity)
rotation = np.zeros((6, 6))
rotation[:3, :3] = rotation[3:, 3:] = [[1., 0., 0.], [0., c, s], [0., -s, c]]
sun_state = np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz])
initial_ecl = rotation @ (initial - sun_state)
command = [args.trajecta, de440, sb441] + [format(x, ".17g") for x in initial_ecl]
raw = subprocess.check_output(command, text=True)
rows = np.array([[float(value) for value in line.split()] for line in raw.splitlines()])
if rows.shape != (3653, 7):
    raise RuntimeError(f"unexpected Trajecta output shape {rows.shape}")
excluded_raw = subprocess.check_output(command + ["--exclude-ceres"], text=True)
excluded = np.array([[float(value) for value in line.split()] for line in excluded_raw.splitlines()])
if excluded.shape != rows.shape or not np.array_equal(excluded[:, 0], rows[:, 0]):
    raise RuntimeError("Ceres-excluded output does not match the reference sample grid")

sim = rebound.Simulation()
sim.add(x=initial[0], y=initial[1], z=initial[2], vx=initial[3], vy=initial[4], vz=initial[5])
extras = assist.Extras(sim, ephem)
extras.gr_eih_sources = 11
extras.forces = [f for f in extras.forces if f not in ("EARTH_HARMONICS", "NON_GRAVITATIONAL")]
positions, velocities, distances_assist, distances_trajecta, earth_distances = [], [], [], [], []
assist_samples = []
au_m = 149597870700.0
for row in rows:
    day = row[0]
    extras.integrate_or_interpolate(day)
    p = sim.particles[0]
    sun = ephem.get_particle("Sun", day)
    ceres = ephem.get_particle("ceres", day)
    earth = ephem.get_particle("Earth", day)
    state = rotation @ (np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz]) -
                        np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz]))
    assist_samples.append(state)
    ceres_ecl = rotation[:3, :3] @ np.array([ceres.x - sun.x, ceres.y - sun.y, ceres.z - sun.z])
    positions.append(float(np.linalg.norm(row[1:4] - state[:3]) * au_m))
    velocities.append(float(np.linalg.norm(row[4:7] - state[3:]) * au_m / 86400.0))
    distances_assist.append(float(np.linalg.norm(np.array([p.x-ceres.x, p.y-ceres.y, p.z-ceres.z]))))
    distances_trajecta.append(float(np.linalg.norm(row[1:4] - ceres_ecl)))
    earth_distances.append(float(np.linalg.norm(np.array([p.x-earth.x, p.y-earth.y, p.z-earth.z]))))

index_assist = int(np.argmin(distances_assist))
index_trajecta = int(np.argmin(distances_trajecta))
index_max = int(np.argmax(positions))
output = {"assist_version": assist.__version__, "rebound_version": rebound.__version__,
          "forces": extras.forces, "gr_eih_sources": extras.gr_eih_sources,
          "epsilon": sim.ri_ias15.epsilon, "min_dt_days": sim.ri_ias15.min_dt,
          "adaptive_mode": str(sim.ri_ias15.adaptive_mode), "sample_count": len(rows),
          "minimum_sampled_earth_distance_au": min(earth_distances),
          "assist_closest_jd_tdb": float(ephem.jd_ref + rows[index_assist, 0]),
          "assist_closest_au": distances_assist[index_assist],
          "trajecta_closest_jd_tdb": float(ephem.jd_ref + rows[index_trajecta, 0]),
          "trajecta_closest_au": distances_trajecta[index_trajecta],
          "difference_at_assist_closest_m": positions[index_assist],
          "max_position_difference_m": positions[index_max],
          "max_difference_jd_tdb": float(ephem.jd_ref + rows[index_max, 0]),
          "final_position_difference_m": positions[-1],
          "final_velocity_difference_m_s": velocities[-1],
          "sampled_ceres_distance_difference_m": abs(distances_assist[index_assist] - distances_trajecta[index_assist]) * au_m,
          "trajecta_ceres_exclusion_effect_m_at_approach": float(np.linalg.norm(rows[index_assist, 1:4] - excluded[index_assist, 1:4]) * au_m),
          "trajecta_ceres_exclusion_effect_m_at_final": float(np.linalg.norm(rows[-1, 1:4] - excluded[-1, 1:4]) * au_m)}
if args.output_schedule_check or args.ias15_check or args.ra15_check:
    def final_state(options):
        run = subprocess.run(command + options, check=True, text=True, capture_output=True)
        values = np.array([float(value) for value in run.stdout.split()])
        if values.shape != (7,):
            raise RuntimeError("unexpected final-state output")
        return values[1:7], run.stderr.strip()

    adaptive_grid = rows[-1, 1:7]
    def difference(left, right):
        delta = left - right
        return {"position_m": float(np.linalg.norm(delta[:3]) * au_m),
                "velocity_m_s": float(np.linalg.norm(delta[3:]) * au_m / 86400.)}

if args.ra15_check:
    reference = rotation @ (np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz]) -
                            np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz]))
    fixed, fixed_report = final_state(["--fixed-step", "1"])
    output["ra15_initial_step_checks"] = []
    output["ra15_final_state_au_au_per_day"] = rows[-1, 1:7].tolist()
    output["ra15_without_ceres_final_state_au_au_per_day"] = excluded[-1, 1:7].tolist()
    output["ra15_fixed_one_day_report"] = fixed_report
    output["ra15_sampled_differences"] = [
        {"requested_day": day, "sample_day": float(rows[index, 0]),
         "vs_assist_m": positions[index]}
        for day in (10, 100, 500, 600, 610, 618, 619, 620, 630, 700, 1000, 2000, 3653)
        for index in (int(np.argmin(np.abs(rows[:, 0] - day))),)]
    for step in (0.001, 0.01, 0.1, 1., 2., 5., 10., 20.):
        state, report = final_state(["--ra15-initial-step", str(step)])
        output["ra15_initial_step_checks"].append({"initial_step_days": step,
            "vs_assist": difference(state, reference), "vs_fixed_one_day": difference(state, fixed),
            "vs_adaptive_initial_10_days": difference(state, adaptive_grid), "report": report})

if args.output_schedule_check:
    adaptive_final, adaptive_final_report = final_state(["--adaptive-final"])
    fixed_final, fixed_final_report = final_state(["--fixed-step", "1"])
    fixed_grid, fixed_grid_report = final_state(["--fixed-step", "1", "--sample-grid"])

    output["output_schedule_check"] = {
        "adaptive_grid_vs_final_only": difference(adaptive_grid, adaptive_final),
        "fixed_grid_vs_final_only": difference(fixed_grid, fixed_final),
        "adaptive_vs_fixed_both_grid": difference(adaptive_grid, fixed_grid),
        "adaptive_vs_fixed_both_final_only": difference(adaptive_final, fixed_final),
        "adaptive_final_report": adaptive_final_report,
        "fixed_final_report": fixed_final_report,
        "fixed_grid_report": fixed_grid_report}
if args.ias15_check:
    steps_run = subprocess.run(command + ["--ias15-steps"], check=True, text=True, capture_output=True)
    accepted_steps = np.array([[float(x) for x in line.split()[1:]] for line in steps_run.stderr.splitlines()
                               if line.startswith("step ")])
    ias_grid, ias_grid_report = final_state(["--ias15-grid"])
    ias_final, ias_final_report = final_state(["--ias15-final"])
    assist_final = rotation @ (np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz]) -
                               np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz]))
    fixed_ias = {}
    output["ias15_check"] = {
        "adaptive_grid_vs_final": difference(ias_grid, ias_final),
        "adaptive_vs_assist": difference(ias_grid, assist_final),
        "adaptive_vs_ra15": difference(ias_grid, adaptive_grid),
        "adaptive_grid_report": ias_grid_report,
        "adaptive_final_report": ias_final_report,
        "accepted_step_days_percentiles": dict(zip(("min", "p10", "median", "p90", "max"),
                                                     map(float, np.percentile(accepted_steps[:, 1], [0, 10, 50, 90, 100])))),
        "shortest_step_end_day": float(accepted_steps[np.argmin(accepted_steps[:, 1]), 0]),
        "fixed_steps": []}
    trace_run = subprocess.run(command + ["--ias15-trace"], check=True, text=True, capture_output=True)
    trace = np.array([[float(value) for value in line.split()] for line in trace_run.stdout.splitlines()])
    if trace.shape != rows.shape or not np.array_equal(trace[:, 0], rows[:, 0]):
        raise RuntimeError("IAS15 trace does not match the reference sample grid")
    ias_differences = np.linalg.norm(trace[:, 1:4] - np.array(assist_samples)[:, :3], axis=1) * au_m
    at_encounter = float(ias_differences[index_assist])
    output["ias15_check"]["position_difference_at_ceres_sample_m"] = at_encounter
    output["ias15_check"]["sampled_thresholds"] = {}
    for threshold in (0.001, 0.01, 0.1, 1.0):
        indices = np.flatnonzero(ias_differences > threshold)
        output["ias15_check"]["sampled_thresholds"][str(threshold)] = (
            {"day": float(rows[indices[0], 0]), "difference_m": float(ias_differences[indices[0]])}
            if len(indices) else None)
    for step in (8., 4., 2., 1.):
        fixed, report = final_state(["--ias15-fixed-step", str(step)])
        fixed_ias[step] = fixed
        output["ias15_check"]["fixed_steps"].append({"days": step,
            "vs_adaptive": difference(fixed, ias_final), "vs_assist": difference(fixed, assist_final),
            "report": report})
    output["ias15_check"]["one_day_vs_two_day"] = difference(fixed_ias[1.], fixed_ias[2.])
    output["ias15_check"]["initial_step_checks"] = []
    for initial_step in (0.001, 0.01, 0.1, 1., 10.):
        initial_run, report = final_state(["--ias15-initial-step", str(initial_step)])
        output["ias15_check"]["initial_step_checks"].append({"initial_step_days": initial_step,
            "vs_assist": difference(initial_run, assist_final),
            "vs_adaptive_initial_10_days": difference(initial_run, ias_final), "report": report})
if args.convergence:
    output["ra15_fixed_step_checks"] = []
    fixed_states = {}
    for step in (8., 4., 2., 1.):
        run = subprocess.run(command + ["--fixed-step", str(step)], check=True, text=True, capture_output=True)
        state = np.array([float(value) for value in run.stdout.split()])
        if state.shape != (7,):
            raise RuntimeError("unexpected RA15 fixed-step output")
        fixed_states[step] = state[1:7]
        output["ra15_fixed_step_checks"].append({"fixed_step_days_after_first": step,
            "position_difference_from_adaptive_m": float(np.linalg.norm(state[1:4] - rows[-1, 1:4]) * au_m),
            "velocity_difference_from_adaptive_m_s": float(np.linalg.norm(state[4:7] - rows[-1, 4:7]) * au_m / 86400.),
            "integration_report": run.stderr.strip()})
    assist_final = rotation @ (np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz]) -
                               np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz]))
    output["ra15_one_day_vs_two_day_position_m"] = float(
        np.linalg.norm(fixed_states[1.][:3] - fixed_states[2.][:3]) * au_m)
    output["ra15_one_day_vs_assist_default_position_m"] = float(
        np.linalg.norm(fixed_states[1.][:3] - assist_final[:3]) * au_m)
    output["assist_ias15_epsilon_checks"] = []
    assist_states = {}
    for epsilon in (1e-9, 1e-10, 1e-11):
        check = rebound.Simulation()
        check.add(x=initial[0], y=initial[1], z=initial[2],
                  vx=initial[3], vy=initial[4], vz=initial[5])
        check_extra = assist.Extras(check, ephem)
        check_extra.gr_eih_sources = 11
        check_extra.forces = extras.forces
        check.ri_ias15.epsilon = epsilon
        for row in rows:
            check_extra.integrate_or_interpolate(row[0])
        q = check.particles[0]
        final = rotation @ (np.array([q.x, q.y, q.z, q.vx, q.vy, q.vz]) -
                            np.array([sun.x, sun.y, sun.z, sun.vx, sun.vy, sun.vz]))
        assist_states[epsilon] = final
        output["assist_ias15_epsilon_checks"].append({"epsilon": epsilon,
            "position_difference_from_default_m": float(np.linalg.norm(final[:3] - assist_final[:3]) * au_m),
            "velocity_difference_from_default_m_s": float(np.linalg.norm(final[3:] - assist_final[3:]) * au_m / 86400.),
            "accepted_steps": check.steps_done})
    output["assist_epsilon_1e10_vs_1e11_position_m"] = float(
        np.linalg.norm(assist_states[1e-10][:3] - assist_states[1e-11][:3]) * au_m)
print(json.dumps(output, indent=2))
