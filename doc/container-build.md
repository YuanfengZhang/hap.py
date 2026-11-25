Building and Running Containers
==============================

This document describes how to build hap.py from source with Docker or
Singularity, test the image, and run a variant comparison. Both images include
Python 3, hap.py, bcftools, samtools, RTG Tools, and Java.

## Requirements

Build on Linux x86-64 using Docker or Git and Singularity. Network access is
required to retrieve the base image and dependencies. Reference genomes and
benchmark datasets are supplied separately when running a comparison.

## Build with Docker

Clone the Python 3 branch:

```bash
git clone --branch py3 https://github.com/YuanfengZhang/hap.py.git
cd hap.py
```

First, build the dependency image:

```bash
docker build -f Dockerfile.dependencies -t hap.py-dependencies:py3 .
```

Then build hap.py using these dependencies:

```bash
docker build -t hap.py-py3 .
```

The dependency image can be reused when rebuilding hap.py after source changes.
The application build installs hap.py and runs the native unit tests.

To use a different dependency image or version label, pass build arguments:

```bash
docker build \
  --build-arg DEPENDENCY_IMAGE=hap.py-dependencies:py3 \
  --build-arg HAPPY_VERSION=v0.3.16-alpha \
  --build-arg SOURCE_REVISION="$(git rev-parse HEAD)" \
  -t hap.py-py3 .
```

Check the installed version:

```bash
docker run --rm hap.py-py3 hap.py --version
```

Running the image without a command displays hap.py help. To run a tool, specify
its name after the image name:

```bash
docker run --rm hap.py-py3 hap.py --help
docker run --rm hap.py-py3 pre.py --help
docker run --rm hap.py-py3 som.py --help
```

## Test the Docker image

The container test creates a small reference sequence, truth and query VCFs,
and a confident-region BED. No additional test data is required.

Run it from the repository directory:

```bash
mkdir -p test-results/tmp
docker run --rm \
  --user "$(id -u):$(id -g)" --workdir /work \
  -e TMPDIR=/work/tmp -v "$PWD/test-results:/work" \
  hap.py-py3 \
  /opt/python/bin/python3.14 -I /opt/happy-source/src/sh/run_container_smoke.py \
  --work-dir /work/engines
```

The test runs three comparison engines with one and two threads. It checks
variant counts, precision, recall, F1 score, output files, and agreement between
thread counts. A successful run prints:

```text
PASS xcmp threads=1
PASS xcmp threads=2
PASS vcfeval threads=1
PASS vcfeval threads=2
PASS scmp-somatic threads=1
PASS scmp-somatic threads=2
```

Inputs, reports, and per-comparison logs are saved under `test-results/engines`.
The `results.json` file lists completed comparisons. The test exits with a
nonzero status if a comparison or check fails.

The directory passed to `--work-dir` must not already exist. To repeat the test,
choose another directory, such as `/work/engines-2`.

## Prepare input files

Place the reference, truth calls, query calls, and confident regions in a local
`data/` directory:

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

The VCFs must be sorted and BGZF-compressed, with TBI or CSI indexes. The
reference FASTA must have an adjacent `.fai` index. Use the same reference
assembly and matching sequence names for the FASTA, VCFs, and BED. BED coordinates
are zero-based and half-open.

For multi-sample VCFs, select the samples to compare using the options described
in [Hap.py usage](happy.md).

## Run with Docker

Mount the input directory read-only and use a separate writable directory for
reports and temporary files:

```bash
mkdir -p results/scratch
docker run --rm \
  --user "$(id -u):$(id -g)" --workdir /out \
  -e TMPDIR=/out/scratch -e RTG_MEM=4G \
  -v "$PWD/data:/data:ro" -v "$PWD/results:/out" \
  hap.py-py3 \
  hap.py /data/truth.vcf.gz /data/query.vcf.gz \
  -r /data/reference.fa -f /data/confident.bed \
  -o /out/benchmark --engine vcfeval --threads 4 \
  --scratch-prefix /out/scratch -V
```

The first VCF contains truth calls; the second contains query calls. The `-f`
option supplies confident regions, and `-o` gives the output filename prefix.
The `-V` option requests an annotated comparison VCF.

This example uses RTG vcfeval. To use hap.py's default comparison engine, replace
`--engine vcfeval` with `--engine xcmp`. The `scmp-somatic` engine compares allele
presence rather than haplotype equivalence.

The `--user` option writes files with the current user's UID and GID. The host
output and scratch directories must be writable by that user. `TMPDIR` and
`--scratch-prefix` place temporary files in the mounted scratch directory.
`RTG_MEM=4G` sets the Java heap for RTG; allow additional memory for Python and
native tools. Adjust memory and thread counts for the dataset and available
resources.

## Run with Singularity

Build directly from source from the repository root (Git and Singularity required):

```bash
mkdir -p .work/singularity/tmp .work/singularity/cache
sudo env SINGULARITY_TMPDIR="$PWD/.work/singularity/tmp" \
  SINGULARITY_CACHEDIR="$PWD/.work/singularity/cache" \
  singularity build hap.py-py3.sif Singularity.linux64.def
```

The definition downloads build dependencies and compiles hap.py on Linux x86-64.
It does not require Docker or a prebuilt hap.py image. `Singularity.def` provides
the same source build.

Run the comparison with the same input and output mounts:

```bash
mkdir -p results/scratch
singularity run --cleanenv \
  --env TMPDIR=/out/scratch --env RTG_MEM=4G \
  --bind "$PWD/data:/data:ro" --bind "$PWD/results:/out" --pwd /out \
  hap.py-py3.sif \
  /data/truth.vcf.gz /data/query.vcf.gz \
  -r /data/reference.fa -f /data/confident.bed \
  -o /out/benchmark --engine vcfeval --threads 4 \
  --scratch-prefix /out/scratch -V
```

`singularity run` invokes hap.py. Use `singularity exec IMAGE.sif TOOL ...`
to run another installed tool.

Run the SIF test, which generates its own input data:

```bash
mkdir -p test-results/tmp
singularity test hap.py-py3.sif
singularity exec --cleanenv --env TMPDIR=/work/tmp \
  --bind "$PWD/test-results:/work" --pwd /work hap.py-py3.sif \
  /opt/python/bin/python3.14 -I /opt/happy-source/src/sh/run_container_smoke.py \
  --work-dir /work/sif-engines
```

A successful run prints the six `PASS` lines listed above. Use a new work directory
when repeating the test.

## Output files

The example above writes the following files to `results/`:

| File | Description |
|------|-------------|
| `benchmark.summary.csv` | Summary counts, precision, recall, and F1 score |
| `benchmark.extended.csv` | Counts and metrics by variant type and region |
| `benchmark.metrics.json.gz` | Metrics in JSON format |
| `benchmark.roc.all.csv.gz` | Metrics at different quality thresholds |
| `benchmark.vcf.gz` | Annotated comparison VCF, requested with `-V` |

Additional ROC tables are produced for individual variant types and subsets.
See [Hap.py usage](happy.md) for report fields, stratification, and further
command-line options.
