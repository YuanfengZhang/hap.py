#!/usr/bin/env bash
# Explicit provisioning; never called during CMake configuration.
set -euo pipefail
SOURCE=$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")
mkdir -p "$SOURCE/.work/tmp" "$SOURCE/.work/pixi-cache"
TMPDIR="$SOURCE/.work/tmp" PIXI_CACHE_DIR="$SOURCE/.work/pixi-cache" proxychains -q pixi install --locked --manifest-path "$SOURCE/pixi.toml"
TMPDIR="$SOURCE/.work/tmp" PIXI_CACHE_DIR="$SOURCE/.work/pixi-cache" proxychains -q pixi install --locked --manifest-path "$SOURCE/external/toolchain/pixi.toml"
