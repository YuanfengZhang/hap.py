#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Download locked archives, rebuild JNI, and assemble the Linux x86_64 RTG6 distribution."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--java-home', type=Path, required=True)
    parser.add_argument('--cmake', type=Path, required=True)
    parser.add_argument('--cc', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--archive-cache', type=Path)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    root: Path = Path(__file__).resolve().parent
    work: Path = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    results: list[dict] = []

    def run(name: str, command: list[str]) -> None:
        with (work / (name + '.log')).open('w') as log:
            process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        results.append({'name': name, 'command': command, 'exit_code': process.returncode})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if process.returncode:
            raise SystemExit(process.returncode)

    archive_dir: Path = work / 'archives'
    archive_dir.mkdir()
    rows: list[dict] = []
    for filename in ('rtg6-dependencies.lock.json', 'jlibdeflate-sources.lock.json'):
        rows.extend(json.loads((root / filename).read_text())['artifacts'])
    for row in rows:
        directory: Path = archive_dir / row['name']
        directory.mkdir()
        destination: Path = directory / 'artifact'
        cached: Path | None = args.archive_cache / row['name'] / 'artifact' if args.archive_cache else None
        if cached and cached.is_file() and hashlib.sha256(cached.read_bytes()).hexdigest() == row['sha256']:
            shutil.copyfile(cached, destination)
        else:
            run('download-' + row['name'], ['proxychains', '-q', 'curl', '-fL', row['url'], '-o', str(destination)])
        if hashlib.sha256(destination.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('Archive checksum mismatch: ' + row['name'])
    run('native-build', [sys.executable, '-I', str(root / 'build-jlibdeflate.py'),
        '--jni-source', str(archive_dir / 'jlibdeflate/artifact'),
        '--libdeflate-source', str(archive_dir / 'libdeflate/artifact'),
        '--java-jar', str(archive_dir / 'jlibdeflate-0.2.0.jar/artifact'),
        '--java-home', str(args.java_home.resolve()), '--cmake', str(args.cmake.resolve()),
        '--cc', str(args.cc.resolve()), '--jobs', str(args.jobs), '--work', str(work / 'native-build')])
    run('rtg-build', [sys.executable, '-I', str(root / 'build-rtg6.py'),
        '--artifact-dir', str(archive_dir), '--java-home', str(args.java_home.resolve()),
        '--native-build', str(work / 'native-build'), '--work', str(work / 'rtg-build')])
    print(work / 'rtg-build/install')


if __name__ == '__main__':
    main()
