#!/usr/bin/env bash
set -euo pipefail
HARNESS_JOBS=${HARNESS_JOBS-4}
if [[ ! "$HARNESS_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    printf 'HARNESS_JOBS must be a positive decimal integer; got <%s>\n' "$HARNESS_JOBS" >&2
    exit 2
fi
export HARNESS_JOBS
printf 'HARNESS_JOBS=%s\n' "$HARNESS_JOBS"
SOURCE=$(realpath "$1")
BUILD=$(realpath -m "$2")
PREFIX=$(realpath -m "$3")
CC=$(realpath "$4")
mkdir -p "$BUILD/tmp"
export CC
export TMPDIR="$BUILD/tmp"
unset CFLAGS CPPFLAGS LDFLAGS CXXFLAGS LIBRARY_PATH LD_LIBRARY_PATH CPATH
cd "$BUILD"
"$SOURCE/Configure" --prefix="$PREFIX" --libdir=lib shared
make -j16
make test
make install_sw
