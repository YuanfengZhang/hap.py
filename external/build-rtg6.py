#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild RTG/HTSJDK6 offline from pinned archives, retaining legacy strict failures and provenance limitations."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import zipfile


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact-dir', type=Path, action='append', required=True)
    parser.add_argument('--native-build', type=Path, required=True)
    parser.add_argument('--java-home', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    root: Path = Path(__file__).resolve().parent
    lock: dict = json.loads((root / 'rtg6-dependencies.lock.json').read_text())
    artifacts: dict[str, Path] = {}
    for row in lock['artifacts']:
        candidates: list[Path] = [candidate for directory in args.artifact_dir for candidate in
            (directory / row['name'] / 'artifact', directory / row['name']) if candidate.is_file()]
        matches: list[Path] = [candidate.resolve() for candidate in candidates if digest(candidate) == row['sha256']]
        if not matches:
            raise ValueError('Missing checksum-matched offline artifact: ' + row['name'])
        artifacts[row['name']] = matches[0]
    work: Path = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    tmp: Path = work / 'tmp'
    tmp.mkdir()
    env: dict[str, str] = dict(os.environ, TMPDIR=str(tmp), JAVA_HOME=str(args.java_home.resolve()))
    java: Path = args.java_home.resolve() / 'bin'
    results: list[dict] = []

    def run(name: str, command: list[str]) -> None:
        with (work / (name + '.log')).open('w') as log:
            process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env, check=False)
        results.append({'name': name, 'command': command, 'exit_code': process.returncode})
        (work / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        if process.returncode:
            raise SystemExit(process.returncode)

    rtg_parent: Path = work / 'rtg-source'
    rtg_parent.mkdir()
    with tarfile.open(artifacts['rtg-3.13.tar.gz']) as archive:
        archive.extractall(rtg_parent, filter='data')
    rtg: Path = next(rtg_parent.iterdir())
    hts: Path = work / 'hts-source'
    hts.mkdir()
    with zipfile.ZipFile(artifacts['htsjdk-6.0.0-sources.jar']) as archive:
        for name in archive.namelist():
            if Path(name).is_absolute() or '..' in Path(name).parts:
                raise ValueError('Unsafe source archive member')
        archive.extractall(hts)
    for name, tree in [('htsjdk-6.0.0-rtg.patch', hts), ('rtg-3.13-modern-java.patch', rtg),
                        ('rtg-3.13-htsjdk6.patch', rtg)]:
        run(name, ['patch', '--batch', '--fuzz=0', '-p1', '-d', str(tree), '-i', str(root / 'patches' / name)])
    output: Path = work / 'install'
    lib: Path = output / 'lib'
    lib.mkdir(parents=True)
    for name, source in artifacts.items():
        if name.endswith('.jar') and not name.endswith('-sources.jar'):
            shutil.copyfile(source, lib / name)
    run('package-jni', [sys.executable, '-I', str(root / 'package-jlibdeflate.py'),
        '--upstream-jar', str(lib / 'jlibdeflate-0.2.0.jar'), '--native-build',
        str(args.native_build.resolve()), '--work', str(work / 'jni-package')])
    shutil.copyfile(work / 'jni-package/jlibdeflate-0.2.0.jar', lib / 'jlibdeflate-0.2.0.jar')
    for name in ('RPlot.jar', 'gzipfix.jar', 'jumble-annotations.jar'):
        shutil.copyfile(rtg / 'lib' / name, lib / name)
    overlay: Path = work / 'overlay'
    classes: Path = work / 'classes'
    tests: Path = work / 'test-classes'
    for directory in (overlay, classes, tests):
        directory.mkdir()
    changed: list[str] = [line[6:] for line in (root / 'patches/htsjdk-6.0.0-rtg.patch').read_text().splitlines() if line.startswith('+++ b/')]
    overlay_jar: Path = lib / 'htsjdk-6.0.0-rtg-overlay.jar'
    cp: str = str(overlay_jar) + ':' + str(lib / '*')
    for name, paths, target, classpath in [
        ('overlay', [hts / name for name in changed], overlay, str(lib / '*')),
        ('main', sorted((rtg / 'src/main/java').rglob('*.java')), classes, cp),
        ('tests', sorted((rtg / 'src/test/java').rglob('*.java')) + [root / 'patches/ModernJexlTest.java'],
         tests, str(classes) + ':' + cp + ':' + str(rtg / 'testLib/spelling.jar')),
    ]:
        argument_file: Path = work / (name + '.args')
        argument_file.write_text('\n'.join(json.dumps(str(path)) for path in paths))
        run(name, [str(java / 'javac'), '-cp', classpath, '-d', str(target), '@' + str(argument_file)])
        if name == 'overlay':
            run('overlay-jar', [str(java / 'jar'), 'cf', str(overlay_jar), '-C', str(overlay), '.'])
    shutil.copytree(rtg / 'src/main/resources', classes, dirs_exist_ok=True)
    for name, value in {'product.name': 'RTG Tools', 'product.version': '3.13-htsjdk6',
                         'build.version': '3.13+htsjdk6', 'build.time': '2026-10-29'}.items():
        (classes / 'com/rtg' / name).write_text(value + '\n')
    test_cp: str = ':'.join(map(str, [tests, rtg / 'src/test/resources', classes, overlay_jar,
                                    lib / '*', rtg / 'testLib/spelling.jar']))
    run('validation', [sys.executable, '-I', str(root / 'validate-htsjdk6.py'), '--java-home',
        str(args.java_home.resolve()), '--classpath', test_cp, '--work', str(work / 'validation')])
    run('rtg-jar', [str(java / 'jar'), 'cf', str(output / 'RTG.jar'), '-C', str(classes), '.'])
    notices: Path = output / 'licenses'
    shutil.copytree(root / 'rtg-licenses', notices)
    shutil.copytree(work / 'jni-package', notices / 'jni-build', ignore=shutil.ignore_patterns('*.jar'))
    shutil.copyfile(artifacts['htsjdk-6.0.0-sources.jar'], notices / 'htsjdk-6.0.0-sources.jar')
    shutil.copyfile(rtg / 'LICENSE.txt', notices / 'RTG-LICENSE.txt')
    shutil.copyfile(rtg / 'installer/README.txt', notices / 'RTG-UPSTREAM-NOTICES.txt')
    shutil.copyfile(root / 'rtg6-dependencies.lock.json', notices / 'build-sources.lock.json')
    for name in ('htsjdk-6.0.0-rtg.patch', 'rtg-3.13-modern-java.patch', 'rtg-3.13-htsjdk6.patch'):
        shutil.copyfile(root / 'patches' / name, notices / name)
    shutil.copytree(hts, notices / 'overlay-source-context')
    manifest: list[dict] = []
    for artifact in [output / 'RTG.jar', *sorted(lib.glob('*.jar'))]:
        embedded: list[str] = []
        with zipfile.ZipFile(artifact) as archive:
            for member in archive.namelist():
                if not member.endswith('/') and any(token in member.lower() for token in ('license', 'notice', 'copying')) and not member.endswith(('.class', '.java')):
                    if Path(member).is_absolute() or '..' in Path(member).parts:
                        raise ValueError('Unsafe notice member')
                    destination: Path = notices / 'jar-notices' / artifact.name / member
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    destination.write_bytes(archive.read(member))
                    embedded.append(member)
        manifest.append({'name': artifact.name, 'sha256': digest(artifact),
                         'upstream_locked_source': next((row for row in lock['artifacts'] if row['name'] == artifact.name), None),
                         'derived_build_receipt': 'jni-build/identity.json' if artifact.name == 'jlibdeflate-0.2.0.jar' else 'build-receipt.json' if artifact.name in ('RTG.jar', 'htsjdk-6.0.0-rtg-overlay.jar') else None,
                         'embedded_notice_members': embedded})
    (notices / 'components.json').write_text(json.dumps({
        'components': manifest,
        'htsjdk6_evidence': 'Full exact-version source jar, original file headers and overlay-source-context.',
        'historical_material_only': ['htsjdk-5.0.2/README.upstream'],
        'patches': {name: digest(notices / name) for name in
                    ('htsjdk-6.0.0-rtg.patch', 'rtg-3.13-modern-java.patch', 'rtg-3.13-htsjdk6.patch')},
        'unresolved_provenance': ['RPlot', 'jumble', 'overlay file-level authorization',
                                  'snappy compiler and Hadoop CRC32C'],
    }, indent=2) + '\n')
    (notices / 'build-receipt.json').write_text(json.dumps({
        'source_lock_sha256': digest(root / 'rtg6-dependencies.lock.json'),
        'source_artifacts': lock['artifacts'],
        'patch_sha256': {name: digest(root / 'patches' / name) for name in
            ('htsjdk-6.0.0-rtg.patch', 'rtg-3.13-modern-java.patch', 'rtg-3.13-htsjdk6.patch')},
        'builder_sha256': digest(Path(__file__)),
        'steps': results,
    }, indent=2) + '\n')
    launcher: str = (root / 'rtg-launcher.sh').read_text().replace('htsjdk-5.0.2-rtg-overlay.jar', overlay_jar.name)
    (output / 'rtg').write_text(launcher)
    (output / 'rtg').chmod(0o755)
    run('version', [str(output / 'rtg'), 'version'])
    (work / 'identity.json').write_text(json.dumps({
        'lock_sha256': digest(root / 'rtg6-dependencies.lock.json'),
        'distributed_jars': {path.name: digest(path) for path in [output / 'RTG.jar', *sorted(lib.glob('*.jar'))]},
        'production_selected': False,
        'unresolved_provenance': ['RPlot', 'jumble', 'overlay file-level authorization',
                                  'snappy compiler and Hadoop CRC32C'],
        'legacy_strict_pass': False,
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
