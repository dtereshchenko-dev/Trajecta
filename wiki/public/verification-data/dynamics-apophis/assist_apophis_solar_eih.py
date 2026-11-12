import argparse
import json
import math
import os
import re

import assist
import numpy as np
import rebound

parser = argparse.ArgumentParser(description="Propagate Apophis with ASSIST and compare with Horizons vectors.")
parser.add_argument("references", nargs="+", help="Saved Horizons vector responses in JSON")
parser.add_argument("--earth-harmonics", choices=("none", "j2", "j234"), default="none")
parser.add_argument("--sampling", choices=("interpolate", "exact"), default="interpolate")
parser.add_argument("--min-dt", type=float, default=0.001, help="Minimum IAS15 step in days")
parser.add_argument("--epsilon", type=float, default=1e-9, help="IAS15 step-control parameter")
args = parser.parse_args()
if not math.isfinite(args.min_dt) or args.min_dt < 0 or not math.isfinite(args.epsilon) or args.epsilon <= 0:
    parser.error("--min-dt must be finite and non-negative; --epsilon must be finite and positive")

def read_reference(path):
    result = json.load(open(path, encoding="utf-8"))["result"]
    lines = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0].splitlines()
    records = []
    for i, line in enumerate(lines):
        if " = A.D." not in line:
            continue
        mjd = float(line.split("=")[0]) - 2400000.5
        values = []
        for row in lines[i + 1:i + 3]:
            values.extend(float(x) for x in re.findall(r"=\s*([+-]?\d+\.\d+E[+-]\d+)", row))
        if len(values) != 6:
            raise ValueError(mjd)
        records.append((mjd, values[:3], values[3:]))
    return records


records = {r[0]: r for path in args.references for r in read_reference(path)}
refs = [records[key] for key in sorted(records)]
if not refs:
    parser.error("No Horizons vector records found")
ephem = assist.Ephem(os.environ['DE440_FILE'], os.environ['SB441_FILE'])
if args.earth_harmonics == "j2":
    ephem.J3E = 0.0
    ephem.J4E = 0.0
obliquity = math.radians(84381.448 / 3600.0)
c, s = math.cos(obliquity), math.sin(obliquity)


def equatorial(v):
    return [v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]]


def ecliptic(v):
    return [v[0], c * v[1] + s * v[2], -s * v[1] + c * v[2]]


def setup(mjd, state):
    ephem.jd_ref = 2400000.5 + mjd
    sim = rebound.Simulation()
    extras = assist.Extras(sim, ephem)
    sim.ri_ias15.epsilon = args.epsilon
    sim.ri_ias15.min_dt = args.min_dt
    extras.gr_eih_sources = 1
    if 'NON_GRAVITATIONAL' not in extras.forces:
        extras.forces = extras.forces + ['NON_GRAVITATIONAL']
    extras.alpha, extras.nk, extras.nm, extras.nn, extras.r0 = 1., 0., 2., 5.093, 1.
    extras.particle_params = np.array([4.999999873689e-13, -2.901766720242e-14, 0.])
    sim.add(x=state[0], y=state[1], z=state[2], vx=state[3], vy=state[4], vz=state[5])
    return sim, extras


start = refs[0][0]
ephem.jd_ref = 2400000.5 + start
sun = ephem.get_particle('Sun', 0.0)
r, v = equatorial(refs[0][1]), equatorial(refs[0][2])
sim, extras = setup(start, [r[0] + sun.x, r[1] + sun.y, r[2] + sun.z,
                            v[0] + sun.vx, v[1] + sun.vy, v[2] + sun.vz])
if args.earth_harmonics == "none":
    extras.forces = [force for force in extras.forces if force != 'EARTH_HARMONICS']
anchor = start
print('SETUP', 'forces', extras.forces, 'epsilon', sim.ri_ias15.epsilon,
      'min_dt', sim.ri_ias15.min_dt, 'adaptive_mode', sim.ri_ias15.adaptive_mode,
      'sampling', args.sampling, 'continuous simulation', flush=True)
for mjd, expected_r, expected_v in refs:
    t = mjd - anchor
    if args.sampling == 'interpolate':
        extras.integrate_or_interpolate(t)
    else:
        sim.integrate(t)
    p = sim.particles[0]
    sun = ephem.get_particle('Sun', t)
    r = ecliptic([p.x - sun.x, p.y - sun.y, p.z - sun.z])
    v = ecliptic([p.vx - sun.vx, p.vy - sun.vy, p.vz - sun.vz])
    dr = [(a - b) * 149597870700. for a, b in zip(r, expected_r)]
    dv = [(a - b) * 149597870700. / 86400. for a, b in zip(v, expected_v)]
    print('RESULT', mjd, 'dr_m', math.hypot(*dr), 'dv_m_s', math.hypot(*dv),
          'xyz_m', *dr, 'vxyz_m_s', *dv, flush=True)
    # Keep the same simulation across output epochs.
print('STEPS', sim.steps_done, flush=True)
