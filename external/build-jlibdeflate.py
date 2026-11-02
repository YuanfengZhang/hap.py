#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build fixed JNI sources offline and verify the actual loaded native library."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jni-source', type=Path, required=True)
    parser.add_argument('--libdeflate-source', type=Path, required=True)
    parser.add_argument('--java-jar', type=Path, required=True)
    parser.add_argument('--java-home', type=Path, required=True)
    parser.add_argument('--cmake', type=Path, required=True)
    parser.add_argument('--cc', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=6)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    root: Path = Path(__file__).resolve().parent
    lock: dict = json.loads((root / 'jlibdeflate-sources.lock.json').read_text())
    archives: dict[str, Path] = {
        'jlibdeflate': args.jni_source.resolve(),
        'libdeflate': args.libdeflate_source.resolve(),
    }
    for row in lock['artifacts']:
        if sha256(archives[row['name']]) != row['sha256']:
            raise ValueError(f'Source checksum mismatch: {row["name"]}')
    work: Path = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    temporary: Path = work / 'tmp'
    temporary.mkdir()
    env: dict[str, str] = dict(os.environ, JAVA_HOME=str(args.java_home.resolve()),
                               TMPDIR=str(temporary))
    trees: dict[str, Path] = {}
    for name, archive in archives.items():
        destination: Path = work / name
        destination.mkdir()
        with tarfile.open(archive) as stream:
            stream.extractall(destination, filter='data')
        members: list[Path] = list(destination.iterdir())
        if len(members) != 1 or not members[0].is_dir():
            raise ValueError(f'Expected one source directory: {name}')
        trees[name] = members[0]
    shutil.copytree(trees['libdeflate'], trees['jlibdeflate'] / 'native/libdeflate',
                    dirs_exist_ok=True)
    notices: Path = work / 'licenses'
    notices.mkdir()
    for name, filename in [('jlibdeflate', 'LICENSE'), ('libdeflate', 'COPYING')]:
        shutil.copyfile(trees[name] / filename, notices / f'{name}-{filename}')
    shutil.copyfile(root / 'jlibdeflate-sources.lock.json', notices / 'sources.lock.json')
    results: list[dict] = []

    def run(name: str, command: list[str]) -> None:
        with (work / f'{name}.log').open('w') as log:
            result = subprocess.run(command, cwd=root.parent, env=env,
                                    stdout=log, stderr=subprocess.STDOUT, check=False)
        results.append({'step': name, 'command': command, 'exit_code': result.returncode})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if result.returncode:
            raise SystemExit(result.returncode)

    cmake: str = str(args.cmake.resolve())
    java: Path = args.java_home.resolve() / 'bin'
    run('compiler', [str(args.cc.resolve()), '--version'])
    run('cmake-version', [cmake, '--version'])
    run('configure', [cmake, '-S', str(trees['jlibdeflate'] / 'native'),
                      '-B', str(work / 'build'), '-DCMAKE_BUILD_TYPE=Release',
                      '-DCMAKE_C_COMPILER=' + str(args.cc.resolve())])
    run('build', [cmake, '--build', str(work / 'build'), '-j', str(args.jobs)])
    library: Path = work / 'build/libjlibdeflate.so'
    classes: Path = work / 'test-classes'
    classes.mkdir()
    jar: Path = args.java_jar.resolve()
    run('test-compile', [str(java / 'javac'), '-cp', str(jar), '-d', str(classes),
                         str(root / 'patches/NativeRoundtrip.java')])
    run('roundtrip', [str(java / 'java'), '--enable-native-access=ALL-UNNAMED',
                      '-Xlog:library=info', '-Djlibdeflate.library.path=' + str(library.parent),
                      '-cp', f'{classes}:{jar}', 'NativeRoundtrip'])
    log: str = (work / 'roundtrip.log').read_text()
    if not any('Loaded library ' in line and str(library) in line for line in log.splitlines()):
        raise RuntimeError('JVM log did not confirm loading rebuilt library')
    identity: dict = {'sources': lock, 'library': str(library),
                      'library_sha256': sha256(library), 'java_jar_sha256': sha256(jar),
                      'test_source_sha256': sha256(root / 'patches/NativeRoundtrip.java'),
                      'license_files': {p.name: sha256(p) for p in sorted(notices.iterdir())},
                      'native_load_verified': True, 'roundtrips': 15,
                      'scope': 'Standalone Linux build; does not reproduce or replace bundled binaries'}
    (work / 'identity.json').write_text(json.dumps(identity, indent=2) + '\n')
    print(f'PASS: {library}; 15 roundtrips; actual native load verified')


if __name__ == '__main__':
    main()
