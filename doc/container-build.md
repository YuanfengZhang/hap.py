# Python 3 hap.py: container build and usage

## Scope

The local image is `hap.py-py3` (`hap.py-py3:latest`). It includes Python,
native tools, RTG, and Java. Only Linux x86_64 has been validated. The image does
not include a whole-genome reference: mount the FASTA, adjacent `.fai`, input
VCFs, indexes, and BED explicitly. The application version remains
`v0.3.16-alpha`; image naming does not change application versioning.

The rebuilt image passed targeted xcmp, vcfeval, and scmp-somatic acceptance
with one and two threads, without network access and as a non-root user. This
does not establish complete historical golden-suite equivalence, acceptance of
all workflow stratifications, or validation on another physical machine.
The image has not been published to a public registry.

## Transfer an existing image

Export on the build machine:

```bash
docker save hap.py-py3 | gzip -1 > hap.py-py3.tar.gz
```

Copy the archive to the target machine and import it:

```bash
gzip -dc hap.py-py3.tar.gz | docker load
docker run --rm --pull never --network none hap.py-py3 hap.py --version
```

The image has no ENTRYPOINT. When overriding its default command, include
`hap.py` or an absolute executable path. With no command, it displays help.

## Run variant benchmarking

Prepare a directory containing matching inputs:

```text
data/
  reference.fa
  reference.fa.fai
  truth.vcf.gz
  truth.vcf.gz.tbi
  query.vcf.gz
  query.vcf.gz.tbi
  confident.bed
```

The FASTA assembly and sequence names must match the VCFs. Use sorted BGZF VCFs
with TBI/CSI indexes; renaming ordinary gzip files does not make them BGZF.
BED coordinates are zero-based, half-open. Select the intended truth/query
samples explicitly when working with multisample inputs.

```bash
mkdir -p results/tmp
docker run --rm --pull never --network none \
  --user "$(id -u):$(id -g)" --workdir /out \
  -e TMPDIR=/out/tmp -e RTG_MEM=4G \
  -v "$PWD/data:/data:ro" -v "$PWD/results:/out" \
  hap.py-py3 \
  hap.py /data/truth.vcf.gz /data/query.vcf.gz \
  -r /data/reference.fa -f /data/confident.bed \
  -o /out/benchmark --engine vcfeval --threads 4 \
  --scratch-prefix /out/tmp
```

`RTG_MEM` limits the Java heap, not total container memory. Allow additional
memory for native tools and Python. Choose threads, memory, and scratch capacity
for your workload. Keep inputs read-only and outputs writable; never write
outputs into the image. `xcmp` is an alternative engine. `scmp-somatic` performs
allele comparisons, not haplotype-equivalence evaluation.

Outputs include `benchmark.summary.csv`, `benchmark.extended.csv`,
`benchmark.vcf.gz`, `benchmark.roc.*.csv.gz`, and `benchmark.metrics.json.gz`.
Recall is truth TP/(TP+FN); precision is query TP/(TP+FP). Different variant
representations can yield different truth and query TP counts. Sample-level
`BD` values must not be interpreted as mutually exclusive record categories.

## Build from source

Use sufficient Docker storage and build resources. Building dependencies requires
network access; runtime benchmarking does not.

```bash
docker build -f Dockerfile.dependencies -t hap.py-dependencies:py3 .

docker build --network none \
  --build-arg DEPENDENCY_IMAGE=hap.py-dependencies:py3 \
  --build-arg HAPPY_VERSION=v0.3.16-alpha \
  -t hap.py-py3 .
```

Pin reused dependency images by immutable digest when reproducibility is required.
The main Dockerfile does not use the historical CentOS/Python 2 build route.
`Dockerfile.ubuntu-with-tests` is an optional historical-suite runner; configure
`HAPPY_IMAGE` when using a different application image.

Dependency integrity is retained through `pixi.lock`,
`external/toolchain/pixi.lock`, `external/bootstrap-sources.lock.json`,
`external/python314-linux64.requirements.lock`, and RTG lock files. Do not remove
source hashes or disable download verification. The current runtime includes
CPython 3.14.10, HTSlib/bcftools/samtools 1.24, OpenSSL 4.1.0, JDK 27, and
RTG 3.13 with an HTSJDK 6 overlay. The two bootstrap Pixi environments and the
separately built Python installation are different prefixes.

Preserve dependency paths inside the image. Copying `/opt/happy` alone is not a
portable installation: it depends on `/opt/python`, `/opt/happy-deps`, and Java.
Docker save/load preserves these paths.

## Targeted validation

Validate installed CLI engines, thread consistency, truth/query counts, and
independently calculated metrics:

```bash
mkdir -p smoke/tmp
docker run --rm --pull never --network none \
  --user "$(id -u):$(id -g)" --workdir /work \
  -e TMPDIR=/work/tmp -v "$PWD/smoke:/work" \
  hap.py-py3 \
  /opt/python/bin/python3.14 -I /opt/happy-source/src/sh/run_container_smoke.py \
  --work-dir /work/engines
```

Permanent regressions are in `src/python/tests/` and `src/sh/run_*_regression.py`.
Run `bin/test_haplotypes` from a writable directory: randomized tests create
relative-path FASTA files. The historical `src/sh/run_tests.sh` suite requires
writable fixtures and the exact hg19 reference. Skipping tests because a
reference is missing must not be reported as passing.

Historical golden outputs still differ from corrected CONF-union,
insertion-boundary, and scmp-distance semantics. Do not rewrite golden outputs
to hide these differences. Small contracts do not replace full real-data
validation. `external/validate-htsjdk6.py` documents precise historical strict
failures and versioned HTSJDK6 validation semantics, not blanket suite success.
Validate your own references, truth sets, and stratifications before large-scale
use. Windows, macOS, and ARM64 builds have not been validated.

## Singularity / SIF

`Singularity.def` provides a general entrypoint; supported platforms depend on
the selected base image. `Singularity.linux64.def` specifically checks Linux
x86_64. Only the latter SIF variant has been built and validated. A SIF does not
provide native binaries for every architecture; Windows/macOS require a Linux
execution environment.

Both definitions use the locally available `hap.py-py3` Docker image. Import the
Docker archive first; no remote image download is required. Example using
Singularity 4.2.2:

```bash
mkdir -p .work/sif/tmp .work/sif/cache
sudo env SINGULARITY_TMPDIR="$PWD/.work/sif/tmp" \
  SINGULARITY_CACHEDIR="$PWD/.work/sif/cache" TMPDIR="$PWD/.work/sif/tmp" \
  singularity build ../hap.py-py3.sif Singularity.linux64.def
```

Building requires access to the local Docker daemon. Some fakeroot setups lose
supplementary Docker socket groups; root builds may be needed. Runtime does not
require root. Base-image digests in the definitions record build provenance;
verify the tag's image identity before building. An existing renamed SIF retains
its original embedded build metadata.

Copy `hap.py-py3.sif` to the target machine; Docker is not required there:

```bash
singularity run --cleanenv hap.py-py3.sif --version
mkdir -p results/tmp
singularity run --cleanenv \
  --env TMPDIR=/out/tmp --env RTG_MEM=4G \
  --bind "$PWD/data:/data:ro" --bind "$PWD/results:/out" --pwd /out \
  hap.py-py3.sif \
  /data/truth.vcf.gz /data/query.vcf.gz \
  -r /data/reference.fa -f /data/confident.bed \
  -o /out/benchmark --engine vcfeval --threads 4 --scratch-prefix /out/tmp
```

`singularity run IMAGE.sif` invokes hap.py, so do not repeat the command name.
Use `singularity exec IMAGE.sif /opt/happy/bin/TOOL` for other tools.
`--cleanenv` isolates ordinary host environment variables, not networking.
Apptainer has similar options, but validation used Singularity.

## Integrity checks and clean build directories

Dependency downloads, pip/Pixi locks, JNI packaging checks, and base-image
digests retain their integrity protections. Each installer invocation creates a
fresh `build-*` directory under `--scratch-path`; source-tree fingerprints are
not used to select reusable caches. Build trees remain available for diagnosis.

The insertion-boundary adapter compares ordered variant keys, FILTER values,
and TRUTH/QUERY GT, BD, BK, BVT, BLT fields against an explicit compressed JSONL
baseline derived from the independently identified legacy VCF. Only the six
reviewed boundary classifications are changed. Fixed input file hashes remain
mandatory; this is not a general-purpose baseline regeneration mechanism.

RTG provisioning uses `provision-rtg6.py` and `build-rtg6.py`. The superseded
HTSJDK 5 build script and its dedicated dependency lock and patch were removed;
shared patches and legal provenance materials remain. CI printing the FASTA
index hash records provenance; it does not compare against an expected hash.

## Non-root and read-only containers

Docker creates `/data` and `/opt/happy-tmp` with mode `1777`, so arbitrary UIDs
can use the default work and temporary directories without root privileges.
A bind mount replaces the image directory permissions: grant the running UID
write access to the host output directory; never recursively chmod input mounts.

A SIF has a read-only root filesystem. Its environment defaults `TMPDIR` to the
current working directory instead of `/opt/happy-tmp`; run from a writable bound
work directory. An explicit `--env TMPDIR=/out/scratch` takes precedence. Create
that directory before execution. This applies to both `singularity run` and
`singularity exec`. Read-only Docker roots likewise require a writable mounted
work/output directory and an explicit writable `TMPDIR` (or a tmpfs there).
