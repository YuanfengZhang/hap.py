#!/usr/bin/env bash
set -euo pipefail
SOURCE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
TYPE=${1:-Debug}
SETUP=${2:-auto}
PREFIX=${3:-"$PWD/install"}
if [[ "$SETUP" != auto ]]; then
    echo "Legacy host-specific setups are retired; supply a CMake toolchain file instead." >&2
    exit 2
fi
[[ "$TYPE" != install ]] || TYPE=Release
if (( $# >= 3 )); then shift 3; else set --; fi
exec cmake -S "$SOURCE" -B "$PWD" -DCMAKE_BUILD_TYPE="$TYPE" -DCMAKE_INSTALL_PREFIX="$PREFIX" "$@"
