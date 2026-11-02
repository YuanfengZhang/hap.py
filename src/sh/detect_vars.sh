#!/usr/bin/env bash
# Source from tests; explicit reference choices are never replaced.
export HG19="${HGREF:-${HG19:-}}"
export HCDIR="${HCDIR:-$PWD/bin}"
if [[ ! -f "$HCDIR/hap.py" ]]; then
    echo "Cannot find hap.py; set HCDIR to the installation bin directory." >&2
    exit 1
fi
export PATH="$HCDIR:$PATH"
export PYTHON="${PYTHON:-python3.14}"
if ! "$PYTHON" -I -c 'import sys; assert sys.version_info[:2] == (3, 14)'; then
    echo "hap.py requires Python 3.14; set PYTHON to its absolute path." >&2
    exit 1
fi
export TMPDIR="${TMPDIR:-$PWD/.work/tmp}"
mkdir -p "$TMPDIR"
HCVERSION=$("$PYTHON" "$HCDIR/hap.py" --version) || exit 1
export HCVERSION
[[ -n "$HCVERSION" ]] || { echo "Cannot determine hap.py version" >&2; exit 1; }
