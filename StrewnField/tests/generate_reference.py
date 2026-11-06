import argparse
import hashlib
import importlib.util
import inspect
import json
import math
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    import sys
    sys.path.insert(0, str(args.source.parent))
    spec = importlib.util.spec_from_file_location("reference", args.source)
    core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)
    rows = []
    for height in core.ALT_GRID:
        temperature, density = core.us76(float(height))
        rows.append([density, temperature, 20.0 * math.exp(-float(height) / 20.0), -7.0])
    rows = np.asarray(rows)
    atmosphere = dict(kind="1d", grid=core.ALT_GRID, rho=rows[:, 0].copy(), temp=rows[:, 1].copy(),
                      u=rows[:, 2].copy(), v=rows[:, 3].copy())
    terrain = dict(lats=np.array([-90., 90.]), lons=np.array([-180., 180.]), h=np.zeros((2, 2)))
    parameters = dict(density=3300., wind_scale=1., cd_scale=1., sigma_abl=8e-8, shape="sphere")
    lat, lon = math.radians(49.8), math.radians(0.8)
    east, north, up = core.enu_basis(lat, lon)
    n = core.WGS84_A / math.sqrt(1 - core.WGS84_E2 * math.sin(lat)**2)
    cases = []
    refined_source = inspect.getsource(core.integrate_fragment)
    replacements = {
        'tol_abs = (0.01, 0.01, 0.01, 1e-4, 1e-4, 1e-4, 1e-9)':
            'tol_abs = (1e-4, 1e-4, 1e-4, 1e-6, 1e-6, 1e-6, 1e-12)',
        '1e-8 * abs(y[c])': '1e-11 * abs(y[c])',
        '), 5.0)': '), 0.05)',
    }
    for old, new in replacements.items():
        assert refined_source.count(old) == 1, old
        refined_source = refined_source.replace(old, new)
    scope = dict(core.__dict__)
    exec(refined_source, scope)
    refined_integrator = scope['integrate_fragment']
    for name, height, mass, velocity, strength, dark in [
        ("dark_50g", 25000., .05, (1200., 0., -1600.), 4e6, True),
        ("entry_first_breakup", 100000., 500., (9000., -1000., -10500.), 4e6, False),
        ("bright_fragment", 30000., 1., (5500., 0., -6500.), 330e6, False),
        ("mass_cutoff", 30000., .005, (10000., 0., -10000.), 330e6, False),
    ]:
        position = np.array([(n+height)*math.cos(lat)*math.cos(lon), (n+height)*math.cos(lat)*math.sin(lon),
                             (n*(1-core.WGS84_E2)+height)*math.sin(lat)])
        vel = velocity[0]*east + velocity[1]*north + velocity[2]*up
        initial = [*position, *vel, mass]
        cfg = dict(core.CONFIG, luminous_speed=1e99 if dark else 3000.)
        core.init_worker(dict(cfg=cfg, atm=atmosphere, terrain=terrain))
        kind, time, state, _ = core.integrate_fragment(initial, dict(parameters, strength=strength), 0.)
        refined_kind, refined_time, refined_state, _ = refined_integrator(initial, dict(parameters, strength=strength), 0.)
        assert kind == refined_kind
        cases.append(dict(name=name, initial=initial, strength_pa=strength, dark=dark,
                          event={"ground": 0, "breakup": 1, "ablated": 2, "timeout": 3}[kind], time_s=time, final=state,
                          refined_time_s=refined_time, refined_final=refined_state))
    output = dict(source_sha256=hashlib.sha256(args.source.read_bytes()).hexdigest(),
                  scope="Single fragments with frozen synthetic horizontal wind; no network or orbit sampling",
                  refined=dict(relative_tolerance=1e-11, position_tolerance_m=1e-4, velocity_tolerance_mps=1e-6,
                               mass_tolerance_kg=1e-12, max_step_s=0.05),
                  atmosphere=rows.tolist(), cases=cases)
    core.CTX['cfg'] = dict(core.CONFIG)
    parent_state = cases[2]['initial'][:6] + [100.0]
    parent = dict(mass=100.0, mult=7.0, gen=0, strength=4e6)
    prm = dict(parameters, strength=4e6, dust=0.45, mu_l=0.225, beta=0.6)
    grouped = core.group_fragments(core.child_masses(55.0, .225, .6, core.CONFIG), core.CONFIG)
    draws = np.arange(len(grouped), dtype=float) / len(grouped)

    class FrozenDraws:
        def __init__(self):
            self.values = iter(draws)

        def uniform(self, low, high):
            return low + (high-low) * next(self.values)

    aux = core.fragment_rhs(parent_state, prm)[1]
    children = core.make_children(parent_state, parent, prm, aux, FrozenDraws())
    output['fragmentation'] = dict(parent=parent_state, groups=grouped, draws=draws.tolist(),
                                   children=[dict(state=c['state'], strength=c['strength'], mult=c['mult']) for c in children])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
