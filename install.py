#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build hap.py against an explicitly provisioned dependency environment."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def check_python(interpreter: str) -> str:
    """Resolve and validate the selected interpreter without invoking a shell."""
    resolved = shutil.which(interpreter)
    if not resolved:
        raise ValueError(f'Python interpreter not found: {interpreter}')
    subprocess.run([resolved, '-I', '-c',
                    'import sys; assert sys.version_info[:2] == (3, 14), '
                    '"hap.py requires Python 3.14"'], check=True)
    return str(Path(resolved).absolute())


def create_build_directory(source: Path, target: Path, scratch: Path) -> Path:
    """Create a fresh build tree; preserve it for diagnostics after failures."""
    source, target, scratch = source.resolve(), target.resolve(), scratch.resolve()
    if scratch == source or scratch == target or target in scratch.parents:
        raise ValueError('Build scratch must be outside the installation target and differ from source')
    scratch.mkdir(parents=True, exist_ok=True)
    return Path(tempfile.mkdtemp(prefix='build-', dir=scratch))


def main() -> None:
    source = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('targetdir', type=Path)
    parser.add_argument('--python-interpreter', default=sys.executable)
    parser.add_argument('--python', choices=['system', 'virtualenv'], default='system')
    parser.add_argument('--python-virtualenv-dir', type=Path)
    parser.add_argument('--python-virtualenv-update', action='store_true')
    parser.add_argument('--dependency-prefix', type=Path, required=True)
    parser.add_argument('--boost-root', type=Path)
    parser.add_argument('--toolchain-prefix', type=Path, required=True)
    parser.add_argument('--scratch-path', type=Path, default=source / '.work' / 'build')
    parser.add_argument('--configuration', choices=['Debug', 'Release', 'RelWithDebInfo'], default='Release')
    parser.add_argument('--build-processes', type=int, default=min(os.cpu_count() or 1, 16))
    parser.add_argument('--version', default='development')
    parser.add_argument('--with-rtgtools', action='store_true')
    parser.add_argument('--rtg-root', type=Path)
    parser.add_argument('--sge-mode', action='store_true')
    parser.add_argument('--sanitizers', action='store_true')
    parser.add_argument('--no-tests', action='store_true')
    args = parser.parse_args()
    interpreter = check_python(args.python_interpreter)
    prefix = args.dependency_prefix.resolve(strict=True)
    toolchain = args.toolchain_prefix.resolve(strict=True)
    cc = toolchain / 'bin' / 'x86_64-conda-linux-gnu-gcc'
    cxx = toolchain / 'bin' / 'x86_64-conda-linux-gnu-g++'
    subprocess.run([str(cxx), '--version'], check=True)
    target = args.targetdir.resolve()
    build = args.scratch_path.resolve()
    if build == source or build == target:
        parser.error('Build directory must differ from source and installation directories')
    if target.exists() and any(target.iterdir()):
        parser.error('Installation target must be empty; use a new versioned target')
    build = create_build_directory(source, target, build)
    temporary = build / 'tmp'
    temporary.mkdir(exist_ok=True)
    env = dict(os.environ, TMPDIR=str(temporary), PATH=f'{prefix / "bin"}{os.pathsep}{os.environ["PATH"]}')
    # Never inherit activated-environment compiler/linker search paths.
    for key in ('CFLAGS', 'CXXFLAGS', 'CPPFLAGS', 'LDFLAGS', 'CPATH',
                'C_INCLUDE_PATH', 'CPLUS_INCLUDE_PATH', 'LIBRARY_PATH', 'LD_LIBRARY_PATH',
                'CMAKE_PREFIX_PATH', 'CMAKE_ARGS', 'CC', 'CXX'):
        env.pop(key, None)
    if args.python == 'virtualenv':
        if not args.python_virtualenv_dir:
            parser.error('--python virtualenv requires --python-virtualenv-dir')
        venv = args.python_virtualenv_dir.resolve()
        if venv.exists() and not args.python_virtualenv_update:
            parser.error('Virtual environment exists; use --python-virtualenv-update')
        subprocess.run([interpreter, '-m', 'venv', str(venv)], check=True, env=env)
        interpreter = str(venv / 'bin' / 'python')
        subprocess.run(['proxychains', '-q', interpreter, '-I', '-m', 'pip', '--isolated', 'install',
                        '--index-url', 'https://pypi.org/simple', '--require-hashes',
                        '-r', str(source / 'external' / 'python314-linux64.requirements.lock')], check=True, env=env)
    # Validate the complete locked runtime closure, including system mode.
    locked = source / 'external' / 'python314-linux64.requirements.lock'
    subprocess.run([interpreter, '-I', '-c',
                    'import importlib.metadata as m, pathlib, sys; '
                    'rows=[x.split()[0].split("==") for x in pathlib.Path(sys.argv[1]).read_text().splitlines() '
                    'if x and not x.startswith("#")]; '
                    'bad=[(n,v,m.version(n)) for n,v in rows if m.version(n)!=v]; '
                    'assert not bad, bad', str(locked)], check=True, env=env)
    # Check runtime requirements before committing to a native build.
    subprocess.run([interpreter, '-I', '-c',
                    'import numpy,pandas,scipy,pysam,psutil,bx.intervals.intersection'], check=True, env=env)
    command = [str(prefix / 'bin' / 'cmake'), '--fresh', '-S', str(source), '-B', str(build),
               f'-DCMAKE_C_COMPILER={cc}', f'-DCMAKE_CXX_COMPILER={cxx}',
               f'-DCMAKE_INSTALL_PREFIX={target}', f'-DCMAKE_BUILD_TYPE={args.configuration}',
               f'-DHAPPY_DEPENDENCY_PREFIX={prefix}', f'-DPython3_EXECUTABLE={interpreter}',
               f'-DHAPLOTYPES_VERSION={args.version}',
               f'-DBUILD_VCFEVAL={"ON" if args.with_rtgtools else "OFF"}',
               '-DRTG_ROOT=',
               f'-DUSE_SGE={"ON" if args.sge_mode else "OFF"}',
               f'-DHAPPY_SANITIZERS={"ON" if args.sanitizers else "OFF"}']
    if args.boost_root:
        command.append(f'-DBoost_ROOT={args.boost_root.resolve(strict=True)}')
    if args.with_rtgtools:
        if not args.rtg_root:
            parser.error('--with-rtgtools requires --rtg-root')
        command.extend(['-DBUILD_VCFEVAL=ON', f'-DRTG_ROOT={args.rtg_root.resolve(strict=True)}'])
    subprocess.run(command, check=True, env=env)
    subprocess.run([command[0], '--build', str(build), '--parallel', str(args.build_processes)], check=True, env=env)
    subprocess.run([command[0], '--install', str(build)], check=True, env=env)
    # Only executable entry points get a shebang, never library modules or a venv.
    for script in (target / 'bin').glob('*.py'):
        text = script.read_text(encoding='utf-8').splitlines(keepends=True)
        if text and text[0].startswith('#!'):
            text.pop(0)
        script.write_text(f'#!{interpreter}\n' + ''.join(text), encoding='utf-8')
        script.chmod(script.stat().st_mode | 0o111)
    if not args.no_tests:
        env.update(PYTHON=interpreter, HCDIR=str(target / 'bin'))
        subprocess.run(['bash', str(source / 'src' / 'sh' / 'run_tests.sh')], cwd=target, check=True, env=env)


if __name__ == '__main__':
    main()
