#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create a provenance-tracked Linux x86_64 JNI jar without retaining unreconstructed platform libraries."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream-jar', type=Path, required=True)
    parser.add_argument('--native-build', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    build: Path = args.native_build.resolve()
    receipt: dict = json.loads((build / 'identity.json').read_text())
    library: Path = Path(receipt['library'])
    lock: dict = json.loads((Path(__file__).resolve().parent / 'jlibdeflate-sources.lock.json').read_text())
    if not (receipt['sources'] == lock):
        raise AssertionError('Native source lock differs')
    if not (receipt['native_load_verified'] and receipt['roundtrips'] == 15):
        raise AssertionError("Validation condition failed: receipt['native_load_verified'] and receipt['roundtrips'] == 15")
    if not (sha256(library) == receipt['library_sha256']):
        raise AssertionError('Native library checksum differs')
    if not (sha256(args.upstream_jar) == receipt['java_jar_sha256']):
        raise AssertionError('Upstream jar checksum differs')
    for name, checksum in receipt['license_files'].items():
        if not (sha256(build / 'licenses' / name) == checksum):
            raise AssertionError(name)
    work: Path = args.work.resolve()
    work.mkdir(parents=True, exist_ok=False)
    target: Path = work / args.upstream_jar.name
    member: str = 'native/linux-x86_64/libjlibdeflate.so'
    removed: list[str] = []
    with zipfile.ZipFile(args.upstream_jar) as source, zipfile.ZipFile(target, 'w') as output:
        for entry in source.infolist():
            if entry.filename.startswith('native/'):
                removed.append(entry.filename)
                continue
            if entry.filename.upper().startswith('META-INF/') and entry.filename.upper().endswith(('.SF', '.RSA', '.DSA', '.EC')):
                raise ValueError('Signed upstream jar requires explicit signature handling')
            output.writestr(entry, source.read(entry.filename))
        info = zipfile.ZipInfo(member, date_time=(2026, 10, 29, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        output.writestr(info, library.read_bytes())
    with zipfile.ZipFile(target) as archive:
        if not ([name for name in archive.namelist() if name.startswith('native/')] == [member]):
            raise AssertionError("Validation condition failed: [name for name in archive.namelist() if name.startswith('native/')] == [member]")
        if not (hashlib.sha256(archive.read(member)).hexdigest() == receipt['library_sha256']):
            raise AssertionError("Validation condition failed: hashlib.sha256(archive.read(member)).hexdigest() == receipt['library_sha256']")
    shutil.copytree(build / 'licenses', work / 'licenses')
    shutil.copyfile(build / 'identity.json', work / 'native-build-receipt.json')
    (work / 'identity.json').write_text(json.dumps({
        'platform': 'linux-x86_64', 'upstream_jar_sha256': sha256(args.upstream_jar),
        'derived_jar_sha256': sha256(target), 'native_member': member,
        'native_sha256': receipt['library_sha256'], 'removed_upstream_native_members': removed,
        'source_lock': lock, 'native_build_receipt_sha256': sha256(build / 'identity.json'),
        'scope': 'Linux x86_64 only; Java classes unchanged; not byte-identical upstream reproduction',
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
