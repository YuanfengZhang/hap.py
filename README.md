Haplotype Comparison Tools
==========================

> [!IMPORTANT]
> This is a customized version of hap.py for [nf-core/variantbenchmarking](https://github.com/nf-core/variantbenchmarking), addressing upstream issues and improving compatibility with Python 3 and current dependencies.
> The original version can be found [here](https://github.com/Illumina/hap.py).

Original author: Peter Krusche <pkrusche@illumina.com>

These tools benchmark variant calls against gold-standard truth datasets using
[HTSlib](https://github.com/samtools/htslib). This customized version retains the
truth-first, query-second CLI and adds Python 3 support, compatibility fixes,
invalid-record quarantine, and complete-record deduplication.

```bash
mkdir -p results/scratch
hap.py truth.vcf.gz query.vcf.gz \
  -r reference.fa -f confident.bed -o results/benchmark \
  --engine vcfeval --threads 4 --scratch-prefix results/scratch -V
```

Use `--engine xcmp` for the default haplotype comparison engine. The `vcfeval`
engine requires RTG and Java; the source-built Docker image includes them.
See [Usage](#usage) for input requirements and current output names.

## Contents

- [Motivation](#motivation)
  - [Complex variant comparison](#complex-variant-comparison)
  - [Variant preprocessing](#variant-preprocessing)
  - [Variant counting](#variant-counting)
- [Usage](#usage)
  - [hap.py](#happy)
  - [pre.py](#prepy)
  - [som.py](#sompy)
- [Workflow integration](#workflow-integration)
- [Installation](#installation)
  - [Docker](#docker)
  - [Source installation](#source-installation)
- [System requirements](#system-requirements)

## Motivation
### Complex variant comparison

A major challenge when comparing VCF files for diploid samples is the handling
of complex variant representations. In a VCF file, we describe two haplotype
sequences by means of REF-ALT pairs and genotypes. These variant calls do not
uniquely represent the haplotype sequences: since alignments are always not unique
even when using a fixed set of gap and substitution scores,
different variant calling methods may produce different variant representations.
While some of these representational differences can be handled using
pre-processing of VCF files (e.g. variant trimming and left-shifting), others
cannot be fixed easily.

In addition to comparing VCF records individually, we produce a graph-based representation
of the VCF alleles, create all possible haplotype sequences, and compare
these by alignment / exact matching. Here is an example where this is needed:

*Variant representation 1 (shown in purple in the image below):*

```
CHROM POS   REF  ALT             GT
chrQ  10    G    GTGTGTGCATGCT   0/1
```

*Variant representation 2 (shown in green in the image below):*

```
CHROM POS   REF  ALT             GT
chrQ  16    G    GCATGCT         0/1
chrQ  19    T    TGTGTG          0/1
```

![](doc/rep_ex.PNG)

Both representations in this example are able to produce the same alt sequences,
but we are not able to match them up with standard VCF tools. In particular,
we can see from this example that the second representation actually may allow us
to create two different sets of alt sequences if they are part of unphased
heterozygous variant calls. When we don't know the phasing
of our variants, the insertions could have occurred on different haplotypes when using
representation 2.

With this tool, we can produce all haplotypes sequences by enumerating paths
through a reference graph. By finding the paths / alt alleles that are
consistent between two VCFs files we can produce accurate benchmarking
numbers for comparing a VCF to a gold standard truth set.
See [doc/spec.md](doc/spec.md) for more information.

An alternative method to compare  complex variant calls is implemented in
[RTG vcfeval](https://github.com/RealTimeGenomics/rtg-tools). It is possible
to use vcfeval with hap.py, and to use hap.py only for pre-processing,
stratification and counting.

The comparison method in vcfeval is more sophisticated than ours and can
resolve some corner cases more accurately.
For whole-genome comparisons, the difference between the two benchmarking
methods is small, but when focusing on difficult subsets of the genome or
when using variant calling methods that produce many complex variant calls,
these corner cases can become relevant. Moreover, when benchmarking against
gold-standard datasets that cover difficult regions of the genome (e.g.
[Platinum Genomes](http://www.illumina.com/platinumgenomes/)), the more complicated
subsets of the genome will be responsible for most of the difference between
methods.

### Variant preprocessing

Another component of hap.py is a variant pre-processing method which
deals with complex variant representations and MNPs. When different callers
may represent variants using a different number of VCF records, we should
attempt to count these in a consistent fashion between methods. One example
is the representation of MNVs as individual SNPs vs. as complex variants.

Consider the following case:

*Complex variant representation*:
```
chrQ  16    GGG    TTT         0/1
```

vs.

*Atomized representation*:
```
chrQ  16    G      T         0/1
chrQ  17    G      T         0/1
chrQ  18    G      T         0/1
```

If this variant is a false-positive, the first representation would naively
contribute a single FP record. A variant caller that outputs the second
representation would instead receive a penalty of three FPs for making
the same variant call. Overall, the difference between the two representations
might show significantly when looking at precision levels or false-positive
rates (since these are relative to the total number of query counts, which
use the same representations), but become important when we need to compare
absolute numbers of false-positives. For this case, hap.py can perform a re-alignment
of REF and ALT alleles on the query VCF, and splits the records into atomic
variant alleles to produce more granular counts using [pre.py](doc/normalisation.md).
Left-shifting and trimming are also supported.

### Variant counting

Hap.py includes a module to produce stratified variant counts. Variant types
are determined using a re-alignment of REF and ALT alleles. This is more reliable
than only using allele lengths. Consider the following complex deletion.

```
chr1    201586350       .       CTCTCTCTCT      CA
```

This complex variant call is equivalent to a deletion, followed by a SNP. Our
quantification code will recognize this variant as a deletion and a SNP, and will
count it in both categories (so a TP call for this variant will contribute a
SNP and an INDEL). This effectively deals with variant calling methods that
prefer to combine local haplotypes in the same variant records
(e.g. Freebayes / Platypus), which would otherwise fall into a hard-to-assess
"COMPLEX" variant call category that varies substantially between
different variant calling methods.

```
chr1    201586350       .       CTCTCTCTC       C
chr1    201586359       .       T               A
```

Another feature of the quantification module in hap.py is stratification into
variant sub-types and into genomic regions. For example, precision and recall
can be computed at the same time for all
[GA4GH stratification regions](https://github.com/ga4gh/benchmarking-tools/tree/master/resources/stratification-bed-files),
and for different INDEL lengths (\<5, 7-15, 16+). Hap.py also calculates
het-hom and Ti/Tv ratios for all subsets of benchmarked variants.
Note that all region matching in hap.py is based on reference coordinates
only. One case where this can lead to counterintuitive results is when considering
hompolymer insertions:

```
Reference:

>chrQ
CAAAAA

VCF:
chrQ    1   C   CA  0/1

BED for homopolymers:
1   6
```

In this example, the variant call given above would not be captured by the bed region for the
homopolymers because it is associated with the reference base just before. To account for this,
the bed intervals need to be expanded to include the padding base just before the regions.

Finally, we produce input data for ROC and precision/recall curves. An
[example](doc/microbench.md) is included.

## Usage

Examples assume an installed toolchain with its `bin` directory on `PATH`.
Run `hap.py --help`, `pre.py --help`, or `som.py --help` for the current options.
Do not invoke the source scripts with Python 2.

Prepare a reference FASTA with its adjacent `.fai`, sorted BGZF-compressed VCFs
with TBI/CSI indexes, and a confident-region BED. Reference assembly and contig
names must match the VCFs. BED coordinates are zero-based and half-open. Select
the intended sample before benchmarking multi-sample inputs. Reference genomes
and benchmark datasets are not bundled with the image.

### hap.py

Compare diploid genotypes and haplotypes; the first VCF is truth, the second query:

```bash
mkdir -p results/scratch
hap.py truth.vcf.gz query.vcf.gz \
  -r reference.fa -f confident.bed -o results/benchmark \
  --engine vcfeval --threads 4 --scratch-prefix results/scratch -V
```

- `-r`: reference FASTA; specify it explicitly for reproducible runs.
- `-f`: confident regions; query calls outside them are classified as UNK.
- `-o`: output prefix, not an output directory.
- `--engine`: `xcmp` (default), `vcfeval`, `scmp-somatic`, or `scmp-distance`.
- `-V`: write the annotated comparison VCF; it is not written by default.
- `--roc QUAL`: select the quality field for threshold tables (QUAL is default).
  Use `--no-roc` to disable these tables. Recall and Precision columns can be
  used to plot a precision-recall curve; Precision is not a false-positive rate.

Typical reports for the command above:

```text
results/benchmark.summary.csv
results/benchmark.extended.csv
results/benchmark.metrics.json.gz
results/benchmark.roc.all.csv.gz
results/benchmark.vcf.gz
```

Additional stratified ROC tables and VCF indexes may also be produced. Counts and
metrics depend on the inputs and selected options; historical example numbers
are not expected results for this version. Recall uses `TRUTH.TP / (TRUTH.TP +
TRUTH.FN)`; Precision uses `QUERY.TP / (QUERY.TP + QUERY.FP)`, excluding UNK.
Zero denominators are reported as zero. See [the evaluation guide](doc/happy.md).

To quarantine supported invalid GT/AD records and evaluate the remaining records:

```bash
hap.py truth.vcf.gz query.vcf.gz \
  -r reference.fa -f confident.bed -o results/benchmark \
  --output-invalid --threads 4 --scratch-prefix results/scratch
```

This writes `results/benchmark.truth.invalid.vcf.gz` and
`results/benchmark.query.invalid.vcf.gz`, each with an adjacent `.reasons.tsv`.
Quarantine scans both inputs before filtering or region selection. Excluding a
record does not mask the corresponding position in the other input and can
change benchmark metrics. Parsing, header, compression, and I/O errors still
fail; an input with no valid records also fails. Existing quarantine outputs are
not overwritten. Without this flag, invalid records are not quarantined.

Final preprocessing removes complete duplicate records by default, including
all sample and annotation fields in the comparison. Different genotypes or
annotations are retained. `--keep-duplicates` disables final deduplication;
with `--bcftools-norm`, it retains the historical early allele-based removal,
so it does not guarantee preservation of every original duplicate.

### pre.py

Normalize one input without running evaluation:

```bash
mkdir -p results/scratch
TMPDIR="$PWD/results/scratch" pre.py query.vcf.gz results/query.pre.vcf.gz \
  -r reference.fa --threads 4 --output-invalid
```

The quarantine file is `results/query.pre.vcf.gz.invalid.vcf.gz`, with an adjacent
`.reasons.tsv`. `--output-invalid` and `--keep-duplicates` have the same meanings
as in hap.py. See [normalisation](doc/normalisation.md) for preprocessing controls.

### som.py

Compare allele presence without genotype or haplotype matching:

```bash
som.py truth.vcf.gz query.vcf.gz \
  -r reference.fa -f confident.bed -o results/somatic
```

The main report is `results/somatic.stats.csv`. Calls outside confident regions
are UNK. See [the somatic comparison guide](doc/sompy.md). The hap.py/pre.py
quarantine and deduplication flags are not som.py options.

## Workflow integration

For [nf-core/variantbenchmarking](https://github.com/nf-core/variantbenchmarking),
changing the host `PATH` does not replace workflow task containers. Configure
`HAPPY_HAPPY` and `HAPPY_PREPY` to use an image built from this source revision,
following the workflow's configuration and container-registry rules.

For example, use an absolute shared SIF path in a custom Nextflow configuration:

```groovy
singularity.enabled = true
singularity.autoMounts = true
docker.enabled = false

process {
    withName: 'HAPPY_(HAPPY|PREPY)' {
        container = '/absolute/path/to/custom-happy.sif'
        containerOptions = '--env TMPDIR=.'
    }
}
```

The SIF path is a placeholder, not a supplied artifact. Set any module-specific
version override to the actual output of `hap.py --version`, rather than copying
a build date from documentation. Load this configuration after site settings
with `-c`, keeping the workflow's required inputs and analysis parameters.

Images, inputs, and work directories must be accessible on every execution node.
For Docker, use an image reference available to each node and compatible with
the workflow's registry settings. Nextflow uses `singularity exec`, not the SIF
runscript. Use a writable scratch directory and retain the engine arguments
required by the analysis. Verify task `.command.run` files and version reports;
use a fresh work directory when validating a changed toolchain. Tool-level tests
do not establish full workflow equivalence.

## Installation

### Docker

Use the published image:

```bash
docker pull ghcr.io/yuanfengzhang/hap.py:v0.3.16-alpha
```

To build from source, build the dependency image before the application image:

```bash
docker build -f Dockerfile.dependencies -t custom-happy-dependencies:build .

docker build --network none \
  --build-arg DEPENDENCY_IMAGE=custom-happy-dependencies:build \
  --build-arg HAPPY_VERSION="$(git rev-parse --short HEAD)" \
  -t custom-happy:build .
```

The image has no ENTRYPOINT: specify `hap.py` after the image name when replacing
the default help command. For example, from a directory containing `data/`:

```bash
mkdir -p results/scratch
docker run --rm --pull never --network none \
  --user "$(id -u):$(id -g)" --workdir /out \
  -e TMPDIR=/out/scratch -e RTG_MEM=4G \
  -v "$PWD/data:/data:ro" -v "$PWD/results:/out" \
  ghcr.io/yuanfengzhang/hap.py:v0.3.16-alpha \
  hap.py /data/truth.vcf.gz /data/query.vcf.gz \
  -r /data/reference.fa -f /data/confident.bed \
  -o /out/benchmark --engine vcfeval --threads 4 \
  --scratch-prefix /out/scratch -V
```

Mount the reference FASTA/FAI, VCFs/indexes, and BED under `data/`. A SIF built
using this repository's definitions has a hap.py runscript, so
`singularity run IMAGE.sif TRUTH QUERY ...` takes hap.py arguments directly;
`singularity exec IMAGE.sif hap.py ...` explicitly selects the executable.
Adapt the definition's base image to the image you built before creating a SIF.

### Source installation

The build uses **CPython 3.14**, an explicitly provisioned dependency prefix, and
private modules in `prefix/lib/hap.py`. `external/bootstrap.sh` provisions the
current dependencies; `install.py --help` documents installation options.
Configuration does not download software or execute Git. See the
[build guide](doc/container-build.md) for dependency provisioning details; adapt
example paths and tags to your environment.

With dependencies already provisioned:

```bash
cmake -S . -B .work/build \
  -DHAPPY_DEPENDENCY_PREFIX=/absolute/path/to/dependencies \
  -DPython3_EXECUTABLE=/absolute/path/to/python3.14 \
  -DHAPLOTYPES_VERSION="$(git rev-parse --short HEAD)" \
  -DCMAKE_BUILD_TYPE=Release
cmake --build .work/build --parallel 4
cmake --install .work/build --prefix "$PWD/install"
```

Preserve the absolute dependency prefixes: moving the installation directory
alone is insufficient. For a sanitizer build, use a separate build directory
with `-DCMAKE_BUILD_TYPE=Debug -DHAPPY_SANITIZERS=ON`. Known historical golden-suite
differences remain; a successful build or targeted test run does not imply that
all historical integration tests pass.

## System requirements

CMake >=3.24, a C++17 compiler, Boost 1.92, HTSlib/bcftools/samtools 1.24,
JsonCpp 1.9.8, zlib 1.3.2, and the packages in `happy.requirements.txt`.
The tested platform is Linux x86-64; historic claims for CentOS 5/6 or Ubuntu
12.04 do not establish support for this toolchain. Windows and macOS have not
been validated during this migration.

A reference need not be human or use chr-prefixed names. The FASTA/FAI sequence
dictionary, not an hg19/hg38 label, is authoritative. Reference precedence is
explicit `-r`, then `HGREF`, then `HG19`, then the documented fallback; an invalid
explicit reference is an error. Full legacy integration tests require the exact
matching hg19 reference, not an arbitrary genome with similar contig names.
