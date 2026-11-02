#!/usr/bin/env bash
set -euo pipefail
ROOT=$(dirname "$(realpath "$0")")
: "${JAVA_HOME:?Set JAVA_HOME to the managed JDK 27 installation}"
: "${TMPDIR:?Set TMPDIR to a writable project-local temporary directory}"
mkdir -p "$TMPDIR"
overlays=("$ROOT"/lib/htsjdk-*-rtg-overlay.jar)
if [[ ${#overlays[@]} -ne 1 || ! -f "${overlays[0]}" ]]; then
    echo "Expected exactly one HTSJDK RTG overlay in $ROOT/lib" >&2
    exit 1
fi
exec "$JAVA_HOME/bin/java" -Xmx"${RTG_MEM:-4g}" -Djava.io.tmpdir="$TMPDIR" -Dtalkback=false -Drtg.version.check-latest=false -cp "$ROOT/RTG.jar:${overlays[0]}:$ROOT/lib/*" com.rtg.RtgTools "$@"
