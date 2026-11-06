import argparse
import copy
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def python_weather_check(snapshot, research):
    import numpy as np
    from PIL import Image

    raw = read(research / 'data/gfs_noaa_20230213_00z_f003_49.5_50.25_-0.5_1.25.json')
    atm = snapshot['atmosphere']
    actual = np.array(atm['samples']).reshape(len(atm['latitudes_deg']), len(atm['longitudes_deg']), -1, 5)
    heights = atm['first_height_m'] + np.arange(atm['height_count']) * atm['height_step_m']
    check = heights <= 40000
    with Image.open(research / 'data/us_nga_egm08_25.tif') as image:
        dx, dy, _ = image.tag_v2[33550]
        i0, j0, _, x0, y0, _ = image.tag_v2[33922]
        keys = image.tag_v2[34735]
        raster_type = next(keys[k+3] for k in range(4, len(keys), 4) if keys[k] == 1025)
        offset = 0 if raster_type == 2 else .5
        geoid = np.asarray(image)
    maximum = np.zeros(5)
    columns = 0
    for a, lat in enumerate(atm['latitudes_deg']):
        if lat not in raw['lats']:
            continue
        i = raw['lats'].index(lat)
        for b, lon in enumerate(atm['longitudes_deg']):
            if lon not in raw['lons']:
                continue
            j = raw['lons'].index(lon)
            x, y = (lon-x0)/dx+i0-offset, (y0-lat)/dy+j0-offset
            ix, iy = int(x), int(y)
            fx, fy = x-ix, y-iy
            n = (1-fy)*((1-fx)*geoid[iy, ix]+fx*geoid[iy, ix+1]) + fy*((1-fx)*geoid[iy+1, ix]+fx*geoid[iy+1, ix+1])
            col = []
            for level in raw['levels']:
                key = str(level)
                h, t = raw['HGT'][key][i][j], raw['TMP'][key][i][j]
                col.append([6356766*h/(6356766-h)+n, level*100/(287.0531*t), t,
                            raw['UGRD'][key][i][j], raw['VGRD'][key][i][j], 0])
            col = np.array(sorted(col))
            expected = np.column_stack([np.exp(np.interp(heights[check], col[:, 0], np.log(col[:, 1])))] +
                                       [np.interp(heights[check], col[:, 0], col[:, k]) for k in range(2, 6)])
            errors = np.abs(actual[a, b, check] - expected)
            errors[:, 0] /= expected[:, 0]
            maximum = np.maximum(maximum, errors.max(axis=0))
            columns += 1
    assert columns >= 20, columns
    assert np.all(maximum < np.array([2e-6, 2e-4, 2e-4, 2e-4, 1e-15])), maximum
    return {'columns': columns, 'heights_per_column': int(check.sum()),
            'max_density_relative_error': maximum[0], 'max_temperature_error_k': maximum[1],
            'max_east_error_mps': maximum[2], 'max_north_error_mps': maximum[3],
            'max_up_error_mps': maximum[4],
            'reference': 'Previously decoded NOAA GFS JSON (rounded to 1e-4), with explicit EGM2008 height shift'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--research', type=Path)
    parser.add_argument('--historical', action='store_true')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[2]
    executable = args.executable.resolve()
    config = read(repo / 'StrewnField/examples/cx1-weather.json')
    source = (repo / 'StrewnField/examples' / config['input']['path']).resolve()
    solution = root / 'work/2026-10-10/2023 CX1'
    solution.mkdir(parents=True)
    shutil.copyfile(source / 'orbit_cartesian_oef.txt', solution / 'orbit_cartesian_oef.txt')
    config['input']['path'] = str(solution)

    def invoke(name, job):
        path = root / f'{name}.json'
        path.write_text(json.dumps(job, indent=2), encoding='utf-8')
        with (root / f'{name}.log').open('w', encoding='utf-8') as log:
            run = subprocess.run([str(executable), str(path)], stdout=log, stderr=subprocess.STDOUT, timeout=900)
        assert run.returncode == 0, (name, run.returncode, (root / f'{name}.log').read_text())
        print(name + ': complete', flush=True)

    invoke('online', config)
    cache = solution / 'StrewnField'
    first = max((cache / 'runs').iterdir())
    summary = read(first / 'summary.json')
    snapshot = read(first / 'environment.json')
    assert summary['completed'] == 200 and summary['workers'] == 15
    assert summary['failed'] == summary['truncated'] == 0
    assert not snapshot['idealized_fallback']
    weather = snapshot['preparation']['weather']
    assert weather['cycle_mjd_utc'] == 59988 and weather['valid_mjd_utc'] == 59988.125
    assert weather['grid_step_deg'] == [.25, .25]
    source_records = [read(path) for folder in ('weather', 'terrain') for path in (cache / folder).glob('*.json')]
    assert all('open-meteo' not in row['source_url'] for row in source_records)
    frozen_hashes = {str(path.relative_to(cache)): hashlib.sha256(path.read_bytes()).hexdigest()
                     for folder in ('weather', 'terrain') for path in (cache / folder).iterdir()}
    config['environment']['offline'] = True
    invoke('offline', config)
    second = max((cache / 'runs').iterdir())
    names = ('entries.json', 'environment.json', 'impacts.csv', 'breakups.jsonl', 'runs.jsonl')
    exact = {name: (first / name).read_bytes() == (second / name).read_bytes() for name in names}
    assert all(exact.values()), exact
    assert all(hashlib.sha256((cache / name).read_bytes()).hexdigest() == checksum for name, checksum in frozen_hashes.items())
    replay = read(first / 'job.json')
    replay['input']['path'] = str(first / 'entries.json')
    replay['environment']['path'] = str(first / 'environment.json')
    replay['output_directory'] = str(root / 'replay')
    invoke('replay', replay)
    assert all((first / name).read_bytes() == (root / 'replay' / name).read_bytes() for name in names)
    report = {'summary': summary, 'run_directory': str(first), 'offline_exact': exact,
              'replay_exact': True, 'source_files': len(source_records),
              'download_bytes': sum(row['bytes'] for row in source_records),
              'gfs_cycle_mjd_utc': weather['cycle_mjd_utc'], 'gfs_valid_mjd_utc': weather['valid_mjd_utc'],
              'gfs_grid_step_deg': weather['grid_step_deg'], 'no_open_meteo_requests': True,
              'orbital_core_sha256': hashlib.sha256((executable.parent / 'OrbitalCore.dll').read_bytes()).hexdigest()}
    previous = repo / 'diagnostics/solution_import_20261010/cli-check/cx1-from-work-final/entries.json'
    report['orbital_entries_byte_exact'] = (first / 'entries.json').read_bytes() == previous.read_bytes()
    assert report['orbital_entries_byte_exact']
    if args.research:
        report['python_weather_comparison'] = python_weather_check(snapshot, args.research)
    if args.historical:
        job = {'schema_version': 1, 'input': {'kind': 'fixed_state', 'state': {
            'epoch_mjd_utc': 58271.75, 'latitude_deg': -21.25, 'longitude_deg': 24.25,
            'height_wgs84_m': 5000, 'velocity_enu_mps': [0, 0, -100]}},
            'physics': {'phase': 'dark_flight'}, 'monte_carlo': {'runs': 1, 'workers': 1},
            'environment': {'mode': 'automatic'}, 'output_directory': str(root / 'ncei-la-date')}
        invoke('ncei-la-date', job)
        old = read(root / 'ncei-la-date/environment.json')['preparation']['weather']
        assert 'oa/prod-model' in old['source_url'] and old['grid_step_deg'] == [.5, .5]
        report['historical_grib2_jpeg2000'] = old
    (root / 'verification.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
