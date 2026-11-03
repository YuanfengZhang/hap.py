# First build Dockerfile.dependencies; DEPENDENCY_IMAGE may be a registry digest.
ARG DEPENDENCY_IMAGE=hap.py-dependencies:py3
FROM ${DEPENDENCY_IMAGE}
ARG HAPPY_VERSION=v0.3.16-alpha
ARG SOURCE_REVISION=uncommitted
LABEL org.opencontainers.image.title="hap.py" \
      org.opencontainers.image.source="https://github.com/YuanfengZhang/hap.py" \
      org.opencontainers.image.version="${HAPPY_VERSION}" \
      org.opencontainers.image.revision="${SOURCE_REVISION}"
ENV JAVA_HOME=/opt/jdk TMPDIR=/opt/happy-tmp
ENV PATH=/opt/happy/bin:/opt/happy-deps/bin:/opt/python/bin:$PATH
WORKDIR /opt/happy-source
COPY . .
COPY external/rtg-licenses/ /opt/rtg/licenses/
RUN mkdir -p /data /opt/happy-tmp && chmod 1777 /data /opt/happy-tmp && \
    for module in /opt/python/lib/python3.14/lib-dynload/_ssl*.so /opt/python/lib/python3.14/lib-dynload/_hashlib*.so; do \
      /opt/happy-toolchain/bin/patchelf --set-rpath /opt/openssl/lib:/opt/happy-deps/lib "$module"; \
    done && \
    /opt/python/bin/python3.14 -I -c 'import ssl,hashlib; print(ssl.OPENSSL_VERSION)' && \
    /opt/python/bin/python3.14 -I install.py /opt/happy \
      --python-interpreter /opt/python/bin/python3.14 \
      --dependency-prefix /opt/happy-deps --toolchain-prefix /opt/happy-toolchain \
      --rtg-root /opt/rtg --with-rtgtools --scratch-path /opt/happy-build \
      --version "$HAPPY_VERSION" --no-tests && \
    mkdir -p /opt/happy-build/test-tmp && \
    bin_dir=/opt/happy/bin && TMPDIR=/opt/happy-build/test-tmp "$bin_dir/test_haplotypes" && \
    "$bin_dir/bcftools" --version && "$bin_dir/hap.py" --version && \
    /opt/happy/libexec/rtg-tools-install/rtg version && \
    rm -rf /opt/happy-build
WORKDIR /data
CMD ["/opt/happy/bin/hap.py", "--help"]
