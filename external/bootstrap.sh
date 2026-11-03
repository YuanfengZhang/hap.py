#!/usr/bin/env bash
# Provision durable, versioned prefixes; never run during CMake configuration.
set -euo pipefail
HARNESS_JOBS=${HARNESS_JOBS-4}
if [[ ! "$HARNESS_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    printf 'HARNESS_JOBS must be a positive decimal integer; got <%s>\n' "$HARNESS_JOBS" >&2
    exit 2
fi
export HARNESS_JOBS
printf 'HARNESS_JOBS=%s\n' "$HARNESS_JOBS"
if [[ $# -lt 1 || $# -gt 2 ]] || [[ $# == 2 && "$2" != --reuse-openssl ]]; then
    printf 'Usage: %s PREFIX [--reuse-openssl]\n' "$0" >&2
    exit 2
fi
ROOT=$(realpath "$(dirname "$0")/..")
PREFIX=$(realpath -m "$1")
ROOT_PREFIX="$PREFIX"
REUSE_OPENSSL=0
if [[ $# == 2 ]]; then
    REUSE_OPENSSL=1
fi
PYTHON_BOOTSTRAP=$(command -v "${PYTHON_BOOTSTRAP:-python3}")
mkdir -p "$PREFIX/tmp"
export TMPDIR="$PREFIX/tmp"
export PIXI_CACHE_DIR="$PREFIX/pixi-cache"
bash "$ROOT/external/make_dependencies.sh"
ARCHIVE_CACHE="$PREFIX/source-bundle"
PREFIX=$(mktemp -d "$PREFIX/run-XXXXXXXX")
printf 'RUN_PREFIX=%s\n' "$PREFIX"
if [[ -e "$ARCHIVE_CACHE" || -L "$ARCHIVE_CACHE" ]]; then
    "$PYTHON_BOOTSTRAP" -I "$ROOT/external/bootstrap-sources.py" "$PREFIX/source-bundle" --archive-cache "$ARCHIVE_CACHE"
else
    "$PYTHON_BOOTSTRAP" -I "$ROOT/external/bootstrap-sources.py" "$PREFIX/source-bundle"
fi
TOOLCHAIN="$ROOT/external/toolchain/.pixi/envs/default"
CC="$TOOLCHAIN/bin/x86_64-conda-linux-gnu-gcc"
export PATH="$TOOLCHAIN/bin:$PATH"
for tool in make patch pkg-config; do
    test -x "$TOOLCHAIN/bin/$tool"
    printf "Build tool: %s\n" "$(command -v "$tool")"
done
OPENSSL_SOURCE=$(echo "$PREFIX"/source-bundle/sources/openssl/*)
PYTHON_SOURCE=$(echo "$PREFIX"/source-bundle/sources/python/*)
JAVA_HOME=$(echo "$PREFIX"/source-bundle/sources/jdk/*)
export JAVA_HOME
SSL_PREFIX="$PREFIX/openssl"
if [[ "$REUSE_OPENSSL" == 1 ]]; then
    "$PYTHON_BOOTSTRAP" -I "$ROOT/external/bootstrap-sources.py" "$ROOT_PREFIX" --retest-openssl
    SSL_PREFIX="$ROOT_PREFIX/openssl"
else
    bash "$ROOT/external/build-openssl.sh" "$OPENSSL_SOURCE" "$PREFIX/build-openssl" "$SSL_PREFIX" "$CC"
fi
bash "$ROOT/external/build-python314.sh" "$PYTHON_SOURCE" "$SSL_PREFIX" "$PREFIX/build-python" "$PREFIX/python" "$CC" "$ROOT/.pixi/envs/default"
"$PREFIX/python/bin/python3.14" -I -m pip --isolated install --index-url https://pypi.org/simple --require-hashes -r "$ROOT/external/python314-linux64.requirements.lock"
"$PYTHON_BOOTSTRAP" -I "$ROOT/external/provision-rtg6.py" --work "$PREFIX/rtg6" --java-home "$JAVA_HOME" \
    --cmake "$ROOT/.pixi/envs/default/bin/cmake" --cc "$CC" --jobs "$HARNESS_JOBS"
printf 'Python: %s\nDependencies: %s\nToolchain: %s\nRTG: %s\nJAVA_HOME: %s\n' "$PREFIX/python/bin/python3.14" "$ROOT/.pixi/envs/default" "$TOOLCHAIN" "$PREFIX/rtg6/rtg-build/install" "$JAVA_HOME"
