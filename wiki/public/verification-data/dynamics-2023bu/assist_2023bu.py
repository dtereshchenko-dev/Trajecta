import json
import math
import os
import re
import sys

import assist
import rebound

mode = sys.argv[3] if len(sys.argv) > 3 else "none"
if mode not in ("none", "j2", "j234"):
    raise ValueError("Earth harmonics mode must be none, j2, or j234")
sampling = sys.argv[4] if len(sys.argv) > 4 else "interpolate"
if sampling not in ("interpolate", "exact"):
    raise ValueError("Sampling must be interpolate or exact")
min_dt = float(sys.argv[5]) if len(sys.argv) > 5 else 0.001
if not math.isfinite(min_dt) or min_dt < 0:
    raise ValueError("Minimum step must be finite and non-negative")

def read_reference(path):
    with open(path, encoding="utf-8-sig") as stream:
        result = json.load(stream)["result"]
    lines = result.split("$$SOE", 1)[1].split("$$EOE", 1)[0].splitlines()
    records = []
    for i, line in enumerate(lines):
        if " = A.D." not in line:
            continue
        mjd = float(line.split("=")[0]) - 2400000.5
        values = []
        for row in lines[i + 1 : i + 3]:
            values.extend(float(x) for x in re.findall(r"=\s*([+-]?\d+\.\d+E[+-]\d+)", row))
        if len(values) != 6:
            raise ValueError(mjd)
        records.append((mjd, values[:3], values[3:]))
    return records


records = {row[0]: row for path in sys.argv[1:3] for row in read_reference(path)}
refs = [records[key] for key in sorted(records)]
if len(refs) < 2:
    raise ValueError("Missing reference vectors")

eps = math.radians(84381.448 / 3600.0)
c, s = math.cos(eps), math.sin(eps)


def equatorial(v):
    return [v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]]


def ecliptic(v):
    return [v[0], c * v[1] + s * v[2], -s * v[1] + c * v[2]]


ephem = assist.Ephem(os.environ["DE440_FILE"], os.environ["SB441_FILE"])
if mode == "j2":
    ephem.J3E = 0.0
    ephem.J4E = 0.0
ephem.jd_ref = 2400000.5 + refs[0][0]
sim = rebound.Simulation()
sim.ri_ias15.min_dt = min_dt
extras = assist.Extras(sim, ephem)
extras.gr_eih_sources = 1
extras.forces = [force for force in extras.forces if force != "NON_GRAVITATIONAL"]
if mode == "none":
    extras.forces = [force for force in extras.forces if force != "EARTH_HARMONICS"]
sun = ephem.get_particle("Sun", 0.0)
r, v = equatorial(refs[0][1]), equatorial(refs[0][2])
sim.add(x=r[0] + sun.x, y=r[1] + sun.y, z=r[2] + sun.z,
        vx=v[0] + sun.vx, vy=v[1] + sun.vy, vz=v[2] + sun.vz)
print("SETUP", "ASSIST", assist.__version__, "REBOUND", rebound.__version__,
      "forces", extras.forces, "eih_sources", extras.gr_eih_sources,
      "epsilon", sim.ri_ias15.epsilon, "min_dt", sim.ri_ias15.min_dt,
      "adaptive_mode", sim.ri_ias15.adaptive_mode, "sampling", sampling, flush=True)

for mjd, expected_r, expected_v in refs:
    t = mjd - refs[0][0]
    if sampling == "interpolate":
        extras.integrate_or_interpolate(t)
    else:
        sim.integrate(t)
    body = sim.particles[0]
    sun = ephem.get_particle("Sun", t)
    r = ecliptic([body.x - sun.x, body.y - sun.y, body.z - sun.z])
    v = ecliptic([body.vx - sun.vx, body.vy - sun.vy, body.vz - sun.vz])
    dr = [(a - b) * 149597870700.0 for a, b in zip(r, expected_r)]
    dv = [(a - b) * 149597870700.0 / 86400.0 for a, b in zip(v, expected_v)]
    print("RESULT", mjd, "dr_m", math.hypot(*dr), "dv_m_s", math.hypot(*dv),
          "xyz_m", *dr, "vxyz_m_s", *dv, flush=True)
print("STEPS", sim.steps_done, flush=True)
