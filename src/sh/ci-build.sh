#!/usr/bin/env bash
set -euo pipefail
: "${PYTHON314:?Absolute managed Python 3.14 interpreter required}"
: "${DEPENDENCY_PREFIX:?Persistent dependency prefix required}"
: "${TOOLCHAIN_PREFIX:?Locked toolchain prefix required}"
: "${RTG_ROOT:?Audited RTG distribution required}"
: "${JAVA_HOME:?Managed JDK required}"
: "${HGREF:?Mounted FASTA required}"
: "${WORKSPACE:?CI workspace required}"
: "${BUILD_NUMBER:?Unique CI build number required}"
test -s "$HGREF"
test -s "$HGREF.fai"
"$PYTHON314" -I -c 'import sys; from pathlib import Path; rows = [line.split("\t") for line in Path(sys.argv[1]).read_text().splitlines()]; sizes = {r[0]: int(r[1]) for r in rows}; required = {"chr1": 249250621, "chr21": 48129895}; bad = {name: sizes.get(name) for name, length in required.items() if sizes.get(name) != length}; sys.exit("Full hg19 fixture reference required; unexpected contig lengths: " + repr(bad) if bad else 0)' "$HGREF.fai"
export TMPDIR="${TMPDIR:-$WORKSPACE/.work/tmp}"
mkdir -p "$TMPDIR"
probe=$(mktemp "$TMPDIR/happy-ci-write.XXXXXXXX")
rm -- "$probe"
printf 'CI temporary filesystem: %s\n' "$TMPDIR"
stat -f -c 'type=%T block_size=%S available_blocks=%a' "$TMPDIR"
df -Pk "$TMPDIR"
# Provenance record only; this prints an identity, not an expected-value check.
sha256sum "$HGREF.fai"
"$PYTHON314" -I "$WORKSPACE/install.py" "$WORKSPACE/.work/install-$BUILD_NUMBER" --python-interpreter "$PYTHON314" --dependency-prefix "$DEPENDENCY_PREFIX" --toolchain-prefix "$TOOLCHAIN_PREFIX" --with-rtgtools --rtg-root "$RTG_ROOT" --scratch-path "$WORKSPACE/.work/build"
