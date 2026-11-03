#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fetch checksum-pinned archives in isolated directories."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
import stat
import shutil
from pathlib import Path, PurePosixPath
import subprocess
import tarfile


def archive_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def verify_tree(archive: Path, destination: Path) -> None:
    """Compare every extracted member without writing or following directory links."""
    if destination.is_symlink() or not destination.is_dir():
        raise ValueError(f'Unsafe source directory: {destination}')
    expected: dict[str, tarfile.TarInfo] = {}
    directories: set[str] = set()
    with tarfile.open(archive) as source:
        for original in source:
            name = PurePosixPath(original.name)
            if name.is_absolute() or '..' in name.parts or not name.parts:
                raise ValueError(f'Unsafe archive path: {original.name}')
            key = str(name)
            if key in expected:
                raise ValueError(f'Duplicate archive member: {key}')
            member = tarfile.data_filter(original, str(destination))
            if member is None or not (member.isfile() or member.isdir() or member.issym() or member.islnk()):
                raise ValueError(f'Unsupported archive member: {key}')
            expected[key] = member
            directories.update(str(parent) for parent in name.parents if str(parent) != '.')
        for key in directories:
            if key in expected and not expected[key].isdir():
                raise ValueError(f'Non-directory archive parent: {key}')
        # extractall chmods each hardlink after linking it: the last member's
        # filtered mode therefore governs every path sharing that inode.
        identities: dict[str, str] = {}
        final_modes: dict[str, int] = {}
        for key, member in expected.items():
            if member.isfile():
                identities[key] = key
            elif member.islnk():
                target = str(PurePosixPath(member.linkname))
                if target not in identities:
                    raise ValueError(f'Hardlink has no preceding regular target: {key}')
                identities[key] = identities[target]
            else:
                continue
            if member.mode is not None:
                final_modes[identities[key]] = member.mode
        wanted = set(expected) | directories
        actual: set[str] = set()
        for folder, folders, files in os.walk(destination, followlinks=False):
            actual.update(str((Path(folder) / name).relative_to(destination)) for name in folders + files)
        if actual != wanted:
            raise ValueError(f'Source tree mismatch {destination}: missing={sorted(wanted - actual)!r}, extra={sorted(actual - wanted)!r}')
        for key in sorted(wanted):
            path = destination / key
            info = path.lstat()
            member = expected.get(key)
            if member is None or member.isdir():
                if not stat.S_ISDIR(info.st_mode):
                    raise ValueError(f'Not a directory: {path}')
            elif member.issym():
                if not stat.S_ISLNK(info.st_mode) or os.readlink(path) != member.linkname:
                    raise ValueError(f'Symlink mismatch: {path}')
            else:
                if not stat.S_ISREG(info.st_mode):
                    raise ValueError(f'Not a regular file: {path}')
                if member.islnk():
                    target = PurePosixPath(member.linkname)
                    if target.is_absolute() or '..' in target.parts or str(target) not in expected:
                        raise ValueError(f'Unsafe hardlink: {path}')
                    target_path = destination / str(target)
                    if target_path.is_symlink() or not target_path.is_file() or not os.path.samefile(path, target_path):
                        raise ValueError(f'Hardlink mismatch: {path}')
                stream = source.extractfile(member)
                if stream is None:
                    raise ValueError(f'Unreadable archive member: {key}')
                with stream, path.open('rb') as current:
                    while True:
                        chunk = stream.read(1024 * 1024)
                        if current.read(len(chunk) if chunk else 1) != chunk:
                            raise ValueError(f'Content mismatch: {path}')
                        if not chunk:
                            break
            mode = final_modes.get(identities[key]) if key in identities else None
            if mode is not None and stat.S_IMODE(info.st_mode) != mode:
                raise ValueError(f'Permission mismatch: {path}: actual={stat.S_IMODE(info.st_mode):o}, expected={mode:o}')
        roots = {PurePosixPath(key).parts[0] for key in wanted}
        if len(roots) != 1 or not (destination / next(iter(roots))).is_dir():
            raise ValueError(f'Expected exactly one source root: {destination}')


def verify_cache(base: Path, rows: list[dict[str, str]]) -> None:
    """Fail closed on any existing cache discrepancy; never repair it."""
    for relative in ('', 'downloads', 'sources'):
        path = base / relative
        if path.is_symlink() or not path.is_dir():
            raise ValueError(f'Unsafe cache directory: {path}')
    names = {row['name'] for row in rows}
    for folder, wanted in ((base, {'downloads', 'sources'}), (base / 'downloads', names), (base / 'sources', names)):
        if {path.name for path in folder.iterdir()} != wanted:
            raise ValueError(f'Unexpected cache entries: {folder}')
    for row in rows:
        download = base / 'downloads' / row['name']
        if download.is_symlink() or not download.is_dir() or {p.name for p in download.iterdir()} != {'archive'}:
            raise ValueError(f'Unsafe download directory: {download}')
        archive = download / 'archive'
        if archive.is_symlink() or not archive.is_file() or archive_digest(archive) != row['sha256']:
            raise ValueError(f'Archive checksum/type mismatch: {archive}')
        verify_tree(archive, base / 'sources' / row['name'])
        print(f"Verified existing source cache: {row['name']}", flush=True)


def reuse_openssl(base: Path) -> None:
    """Validate in-place artifacts and retest them; never install or relocate."""
    # Trusted local build only: configdata Perl and the Makefile execute code.
    print('Trust boundary: executing our trusted local configdata/Makefile; not a sandbox.', flush=True)
    root = Path(__file__).resolve().parent.parent
    base = base.resolve(strict=True)
    source = base / 'source-bundle'
    row = next(row for row in json.loads((root / 'external/bootstrap-sources.lock.json').read_text()) if row['name'] == 'openssl')
    archive = source / 'downloads/openssl/archive'
    if archive.is_symlink() or archive_digest(archive) != row['sha256']:
        raise ValueError('OpenSSL archive identity mismatch')
    verify_tree(archive, source / 'sources/openssl')
    build = base / 'build-openssl'
    prefix = base / 'openssl'
    toolchain = root / 'external/toolchain'
    cc = toolchain / '.pixi/envs/default/bin/x86_64-conda-linux-gnu-gcc'
    config = json.loads(subprocess.check_output([
        '/usr/bin/perl', '-I', str(build), '-Mconfigdata', '-MJSON::PP',
        '-e', r'print JSON::PP->new->encode(\%configdata::config)'], text=True))
    expected_args = [f'--prefix={prefix}', '--libdir=lib', 'shared']
    if (config['CC'] != str(cc) or config['prefix'] != str(prefix)
            or config['version'] != '4.1.0' or config['libdir'] != 'lib'
            or config['perlargv'] != expected_args
            or (build / config['sourcedir']).resolve() != (source / 'sources/openssl/openssl-4.1.0').resolve()):
        raise ValueError('OpenSSL configdata identity mismatch')
    print('Current toolchain lock SHA256:', archive_digest(toolchain / 'pixi.lock'), flush=True)
    print('Current compiler SHA256:', archive_digest(cc), flush=True)
    subprocess.run([str(cc), '--version'], check=True)
    libraries = ['libssl.so.4', 'libcrypto.so.4']
    paths = [folder / name for folder in (build, prefix / 'lib') for name in libraries]
    before = {path: archive_digest(path) for path in paths}
    for name in libraries:
        if before[build / name] != before[prefix / 'lib' / name]:
            raise ValueError(f'Installed/build library mismatch: {name}')
    # Only this CLI inspection gets a library search override. Python must use rpath.
    inspect_env = dict(os.environ, LD_LIBRARY_PATH=str(prefix / 'lib'))
    version = subprocess.check_output([str(prefix / 'bin/openssl'), 'version', '-a'], env=inspect_env, text=True)
    if not version.startswith('OpenSSL 4.1.0 '):
        raise ValueError('Installed OpenSSL version mismatch')
    print(version, flush=True)
    linked = subprocess.check_output(['/usr/bin/ldd', str(prefix / 'bin/openssl')], env=inspect_env, text=True)
    for name in libraries:
        if f'{name} => {prefix / "lib" / name} ' not in linked:
            raise ValueError(f'Unexpected installed library loading: {linked}')
    print(linked, flush=True)
    print('Retesting current artifacts; historical toolchain byte identity is not proven.', flush=True)
    result = subprocess.run([str(toolchain / '.pixi/envs/default/bin/make'), 'test'], cwd=build)
    changed = [str(path) for path in paths if archive_digest(path) != before[path]]
    if changed:
        raise ValueError(f'make test changed critical libraries: {changed}')
    result.check_returncode()
    print('Reused unchanged OpenSSL installation:', prefix, flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--archive-cache', type=Path)
    parser.add_argument('--retest-openssl', action='store_true')
    args = parser.parse_args()
    if args.retest_openssl:
        reuse_openssl(args.directory)
        return
    root = Path(__file__).resolve().parent
    base = args.directory.absolute()
    rows = json.loads((root / 'bootstrap-sources.lock.json').read_text())
    if base.exists() or base.is_symlink():
        verify_cache(base, rows)
        return
    base.mkdir(parents=True, exist_ok=False)
    for row in rows:
        download = base / 'downloads' / row['name']
        download.mkdir(parents=True)
        archive = download / 'archive'
        if args.archive_cache is not None:
            cached = args.archive_cache / 'downloads' / row['name'] / 'archive'
            if cached.is_symlink() or not cached.is_file() or archive_digest(cached) != row['sha256']:
                raise ValueError(f'Invalid reusable archive: {cached}')
            shutil.copyfile(cached, archive)
        else:
            subprocess.run(['curl', '-fL', row['url'], '-o', str(archive)], check=True)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != row['sha256']:
            raise ValueError('Archive checksum mismatch')
        destination = base / 'sources' / row['name']
        destination.mkdir(parents=True)
        with tarfile.open(archive) as source:
            source.extractall(destination, filter='data')
        verify_tree(archive, destination)
        print(row['name'], next(destination.iterdir()), flush=True)


if __name__ == '__main__':
    main()
