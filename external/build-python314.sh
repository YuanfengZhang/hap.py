#!/usr/bin/env bash
# Build verified CPython/OpenSSL sources outside their extraction directories.
set -euo pipefail
if [[ $# != 6 ]]; then
    echo "Usage: $0 PYTHON_SOURCE OPENSSL_PREFIX BUILD_DIR INSTALL_PREFIX COMPILER DEPENDENCY_PREFIX" >&2
    exit 2
fi
SOURCE=$(realpath "$1")
SSL=$(realpath "$2")
BUILD=$(realpath -m "$3")
PREFIX=$(realpath -m "$4")
CC=$(realpath "$5")
DEPS=$(realpath "$6")
ROOT=$(realpath "$(dirname "$0")/..")
test -f "$DEPS/include/zlib.h"
test -f "$DEPS/lib/libz.so"
ZLIB_VERSION=$(sed -n 's/^#define ZLIB_VERSION "\(.*\)"/\1/p' "$DEPS/include/zlib.h")
[[ "$ZLIB_VERSION" == 1.3.2 ]]
export ZLIB_CFLAGS="-I$DEPS/include"
export ZLIB_LIBS="$DEPS/lib/libz.so -Wl,-rpath,$DEPS/lib"
printf 'ZLIB_PREFIX=%s ZLIB_VERSION=%s\n' "$DEPS" "$ZLIB_VERSION"
mkdir -p "$BUILD/tmp"
export TMPDIR="$BUILD/tmp"
unset CFLAGS CXXFLAGS CPPFLAGS LDFLAGS LIBRARY_PATH LD_LIBRARY_PATH CPATH C_INCLUDE_PATH CPLUS_INCLUDE_PATH
export CC
patch --forward -d "$SOURCE" -p1 < "$ROOT/external/patches/cpython-3.14-openssl4.patch"
patch --forward -d "$SOURCE" -p1 < "$ROOT/external/patches/cpython-openssl4-version-test.patch"
cd "$BUILD"
FLAGS="-I$DEPS/include -DOPENSSL_NO_SSL3 -DOPENSSL_NO_SSL3_METHOD -DOPENSSL_NO_TLS1_METHOD -DOPENSSL_NO_TLS1_1_METHOD -DOPENSSL_NO_TLS1_2_METHOD"
LINK_FLAGS="-L$SSL/lib -Wl,-rpath,$SSL/lib -Wl,-rpath,$DEPS/lib"
CPPFLAGS="$FLAGS" LDFLAGS="$LINK_FLAGS" "$SOURCE/configure" --prefix="$PREFIX" --with-openssl="$SSL" --with-openssl-rpath=auto --with-ensurepip=install
rg -q '^MODULE_ZLIB_STATE=yes$' Makefile
make -j16 CFLAGS= LDFLAGS="$LINK_FLAGS" CPPFLAGS="$FLAGS"
for module in _ssl _hashlib; do
    artifacts=(Modules/"$module".cpython-*.so)
    [[ ${#artifacts[@]} == 1 && -f "${artifacts[0]}" ]]
    dynamic=$(readelf -d "${artifacts[0]}")
    if [[ "$dynamic" == *'[libssl.so.3]'* || "$dynamic" == *'[libcrypto.so.3]'* ]]; then
        printf 'Wrong OpenSSL DT_NEEDED for %s\n%s\n' "$module" "$dynamic" >&2
        exit 1
    fi
    [[ "$dynamic" == *'[libcrypto.so.4]'* ]]
    if [[ "$module" == _ssl ]]; then
        [[ "$dynamic" == *'[libssl.so.4]'* ]]
    fi
    printf 'Verified DT_NEEDED for %s\n%s\n' "$module" "$dynamic"
done
./python -I - "$DEPS" "$ZLIB_VERSION" "$SSL" <<'PY_CHECK'
import pathlib
import sys
import zlib
import _ssl
import _hashlib
prefix = pathlib.Path(sys.argv[1]).resolve()
assert zlib.ZLIB_VERSION == zlib.ZLIB_RUNTIME_VERSION == sys.argv[2]
def mapped_paths(text: str) -> set[pathlib.Path]:
    paths: set[pathlib.Path] = set()
    for line in text.splitlines():
        fields = line.split(None, 5)
        if len(fields) == 6 and fields[5].startswith('/'):
            # Linux maps escapes newline as \012, but leaves spaces and
            # literal backslashes untouched. Do not decode other octal text.
            paths.add(pathlib.Path(fields[5].replace('\\012', '\n')).resolve())
    return paths

mapped = mapped_paths(pathlib.Path('/proc/self/maps').read_text())
loaded = {path for path in mapped if path.name.startswith('libz.so')}
assert loaded and all(path.is_relative_to(prefix / 'lib') for path in loaded), loaded
ssl_prefix = pathlib.Path(sys.argv[3]).resolve()
for name in ('libssl.so.4', 'libcrypto.so.4'):
    paths = {path for path in mapped if path.name == name}
    assert paths == {(ssl_prefix / 'lib' / name).resolve()}, paths
print('Verified OpenSSL loading:', _ssl.OPENSSL_VERSION, ssl_prefix)
assert zlib.decompress(zlib.compress(b'bootstrap')) == b'bootstrap'
print('Verified managed zlib:', zlib.ZLIB_RUNTIME_VERSION, sorted(map(str, loaded)))
PY_CHECK
./python -I -m test test_ssl test_hashlib -v
make install CFLAGS= LDFLAGS="$LINK_FLAGS" CPPFLAGS="$FLAGS"
