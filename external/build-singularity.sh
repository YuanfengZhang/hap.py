#!/usr/bin/env bash
# Build the same source toolchain used by the Docker images.
set -euo pipefail
[[ $(uname -s) == Linux && $(uname -m) == x86_64 ]]
export DEBIAN_FRONTEND=noninteractive TMPDIR=/opt/happy-tmp
mkdir -p "$TMPDIR" /data
chmod 1777 "$TMPDIR" /data
apt-get -o APT::Update::Error-Mode=any update
apt-get install -y --no-install-recommends ca-certificates curl bzip2 xz-utils tar bash coreutils build-essential binutils ripgrep git patch
mkdir -p /opt/pixi-download
curl -fL --retry 3 https://github.com/prefix-dev/pixi/releases/download/v0.70.2/pixi-x86_64-unknown-linux-musl.tar.gz -o /opt/pixi-download/pixi.tar.gz
tar -xzf /opt/pixi-download/pixi.tar.gz -C /usr/local/bin pixi
export PATH=/usr/local/bin:$PATH
cd /opt/happy-source
bash external/make_dependencies.sh
PYTHON_BOOTSTRAP=/opt/happy-source/.pixi/envs/default/bin/python \
    HARNESS_JOBS="${HARNESS_JOBS:-4}" bash external/bootstrap.sh /opt/bootstrap
runs=(/opt/bootstrap/run-*)
[[ ${#runs[@]} == 1 ]]
run=${runs[0]}
chmod 0755 "$run"
ln -s /opt/happy-source/.pixi/envs/default /opt/happy-deps
ln -s /opt/happy-source/external/toolchain/.pixi/envs/default /opt/happy-toolchain
ln -s "$run/python" /opt/python
ln -s "$run/openssl" /opt/openssl
ln -s "$run/rtg6/rtg-build/install" /opt/rtg
jdks=("$run"/source-bundle/sources/jdk/*)
[[ ${#jdks[@]} == 1 ]]
ln -s "${jdks[0]}" /opt/jdk
export JAVA_HOME=/opt/jdk CMAKE_GENERATOR=Ninja
export PATH=/opt/happy/bin:/opt/happy-deps/bin:/opt/python/bin:$PATH
for module in /opt/python/lib/python3.14/lib-dynload/_ssl*.so /opt/python/lib/python3.14/lib-dynload/_hashlib*.so; do
    /opt/happy-toolchain/bin/patchelf --set-rpath /opt/openssl/lib:/opt/happy-deps/lib "$module"
done
mkdir -p /opt/rtg/licenses
cp -a external/rtg-licenses/. /opt/rtg/licenses/
/opt/python/bin/python3.14 -I install.py /opt/happy \
    --python-interpreter /opt/python/bin/python3.14 \
    --dependency-prefix /opt/happy-deps --toolchain-prefix /opt/happy-toolchain \
    --rtg-root /opt/rtg --with-rtgtools --scratch-path /opt/happy-build \
    --version v0.3.16-alpha --no-tests
mkdir -p /opt/happy-build/test-tmp
(cd /opt/happy-build/test-tmp && /opt/happy/bin/test_haplotypes)
/opt/python/bin/python3.14 -I -c 'import ssl,numpy,pandas,pysam; print(ssl.OPENSSL_VERSION)'
/opt/happy/bin/hap.py --version
/opt/happy/libexec/rtg-tools-install/rtg version
rm -rf /opt/happy-build /opt/pixi-download /var/lib/apt/lists/*
