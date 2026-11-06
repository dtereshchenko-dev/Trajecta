import argparse
import copy
import hashlib
import json
import shutil
import subprocess
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--online', action='store_true')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[2]
    original = json.loads((repo / 'StrewnField/examples/cx1-work.json').read_text())
    solution = root / 'work/2026-10-10/2023 CX1'
    solution.mkdir(parents=True)
    source = repo / 'diagnostics/entry_mc_20261010/complete/runtime/work/2026-10-10/Sar2667/orbit_cartesian_oef.txt'
    shutil.copyfile(source, solution / source.name)
    original['input']['path'] = str(solution)
    original.pop('output_directory', None)
    base_env = json.loads((repo / 'diagnostics/solution_import_20261010/cli-check/cx1-from-work-final/environment.json').read_text())
    docs = {
        '/weather': {'schema_version': 1, 'height_datum': 'WGS84', 'atmosphere': base_env['atmosphere'], 'source': 'test atmosphere'},
        '/terrain': {'schema_version': 1, 'height_datum': 'WGS84', 'terrain': base_env['terrain'], 'source': 'test terrain'},
        '/invalid': {'schema_version': 1, 'height_datum': 'MSL', 'terrain': base_env['terrain']},
    }
    requests = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            try:
                self.respond()
            except ConnectionError:
                pass

        def respond(self):
            requests.append(self.path)
            if self.path == '/redirect':
                self.send_response(302)
                self.send_header('Location', '/weather')
                self.end_headers()
                return
            if self.path == '/slow':
                time.sleep(2)
            if self.path in ('/oversize', '/oversize-stream'):
                data = b'x' * 2000000
            elif self.path == '/broken':
                data = b'{"schema_version": 1'
            elif self.path in docs:
                data = json.dumps(docs[self.path]).encode()
            else:
                self.send_error(404)
                return
            self.send_response(200)
            if self.path != '/oversize-stream':
                self.send_header('Content-Length', str(len(data) + (100 if self.path == '/broken' else 0)))
            self.end_headers()
            try:
                self.wfile.write(data)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def log_message(self, *_):
            pass

    def invoke(name, config, expected=0):
        path = root / (name + '.json')
        path.write_text(json.dumps(config, indent=2) + '\n')
        result = subprocess.run([str(args.executable.resolve()), str(path)],
                                capture_output=True, text=True, timeout=120)
        (root / (name + '.log')).write_text(result.stdout + result.stderr, encoding='utf-8')
        assert result.returncode == expected, (name, result.returncode, result.stdout, result.stderr)
        return result.stderr

    def find_runs():
        return sorted((solution / 'StrewnField/runs').iterdir())

    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    worker = threading.Thread(target=server.serve_forever, daemon=True)
    worker.start()
    url = f'http://127.0.0.1:{server.server_port}'
    remote = copy.deepcopy(original)
    remote['environment'] = {'mode': 'remote', 'weather_url': url + '/redirect', 'terrain_url': url + '/terrain'}
    try:
        invoke('remote', remote)
        first = find_runs()[0]
        assert requests == ['/redirect', '/weather', '/terrain'], requests
        assert len(list((solution / 'StrewnField/weather').glob('*.json'))) == 1
        assert len(list((solution / 'StrewnField/terrain').glob('*.json'))) == 1
        cache = next((solution / 'StrewnField/weather').glob('*.json'))
        assert json.loads(cache.read_text())['effective_url'] == url + '/weather'
        summary = json.loads((first / 'summary.json').read_text())
        assert summary['completed'] == 200 and summary['failed'] == summary['truncated'] == 0
        previous = repo / 'diagnostics/solution_import_20261010/cli-check/cx1-from-work-final'
        exact = {name: (first / name).read_bytes() == (previous / name).read_bytes()
                 for name in ('impacts.csv', 'entries.json', 'runs.jsonl', 'breakups.jsonl')}
        assert all(exact.values()), exact

        failures = {}
        for name, route, change, message in [
            ('http-404', '/missing', {}, '404'),
            ('content-limit', '/oversize', {'max_download_bytes': 1000000}, 'max_download_bytes'),
            ('stream-limit', '/oversize-stream', {'max_download_bytes': 1000000}, 'max_download_bytes'),
            ('incomplete', '/broken', {}, 'Download failed'),
            ('timeout', '/slow', {'timeout_seconds': 1}, 'Download failed'),
            ('invalid-schema', '/invalid', {}, 'WGS84'),
        ]:
            config = copy.deepcopy(remote)
            config['input']['clones'] = 1
            config['monte_carlo']['runs'] = 1
            config['environment'].update({'terrain_url': url + route, **change})
            before = set(find_runs())
            error = invoke(name, config, 1)
            assert message in error, error
            new = (set(find_runs()) - before).pop()
            assert not (new / 'impacts.csv').exists()
            assert len(list((solution / 'StrewnField/terrain').glob('*.json'))) == 1
            assert not list((solution / 'StrewnField/terrain').glob('*.part'))
            assert not list((solution / 'StrewnField/terrain').glob('*.lock'))
            failures[name] = True
    finally:
        server.shutdown()
        server.server_close()
        worker.join()

    count = len(requests)
    before = set(find_runs())
    invoke('offline-cache', remote)
    second = (set(find_runs()) - before).pop()
    assert len(requests) == count
    for name in exact:
        assert (first / name).read_bytes() == (second / name).read_bytes(), name
    assert first != second
    replay = subprocess.run([str(args.executable.resolve()), str(first / 'job.json')], capture_output=True, text=True, timeout=120)
    assert replay.returncode == 0, replay.stderr
    assert (first / 'impacts.csv').read_bytes() == (first / 'replay/impacts.csv').read_bytes()

    online = None
    if args.online:
        config = copy.deepcopy(remote)
        config['input']['clones'] = 1
        config['environment']['terrain_url'] = 'https://curl.se/libcurl/c/CURLOPT_NOSIGNAL.html'
        error = invoke('https-transport', config, 1)
        assert 'parse_error' in error and 'Download failed' not in error, error
        online = 'HTTPS transfer succeeded; HTML was correctly rejected as non-JSON'
    report = {'remote_pipeline': summary, 'prior_results_byte_exact': exact,
              'default_output_root': str(solution / 'StrewnField'),
              'offline_cache_exact': True, 'self_contained_replay_exact': True,
              'http_failures_rejected': failures, 'https': online,
              'orbital_core_sha256': hashlib.sha256((args.executable.parent / 'OrbitalCore.dll').read_bytes()).hexdigest()}
    (root / 'verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
