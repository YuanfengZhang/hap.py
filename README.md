Haplotype Comparison Tools
==========================

Peter Krusche <pkrusche@illumina.com>

This is a set of programs based on [htslib](https://github.com/samtools/htslib)
to benchmark variant calls against gold standard truth datasets.

To compare a VCF against a gold standard dataset, use the following commmand line
to perform genotype-level haplotype comparison.

```bash
hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa
```

We also have a script to perform comparisons only based on chromosome, position,
and allele identity. This comparison will not resolve haplotypes and only verify
that the same alleles were observed at the same positions (e.g. for comparison
of somatic callsets).

```bash
som.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa
```

More information can be found below in the [usage section](#usage).

## Contents

* [Motivation](#motivation)
  * [Complex variant comparison](#complex-variant-comparison)
  * [Variant preprocessing](#variant-preprocessing)
  * [Variant counting](#variant-counting)
* [Usage](#usage)
  * [hap.py](#happy)
  * [som.py](#sompy)
* [Using Python 3 hap.py in variantbenchmarking](#using-python-3-happy-in-variantbenchmarking)
* [Installation](#installation)
  * [Docker](#docker-recommended)
  * [Helper script](#helper-script)
  * [Compiling from source with CMake](#compiling-from-source-with-cmake)
* [System requirements](#system-requirements)
  * [Hardware](#hardware)
  * [Linux](#linux)
  * [OS X](#os-x)
  * [Windows](#windows)
  * [Other requirements](#other-requirements)

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

The main two tools are hap.py (diploid precision/recall evaluation) and som.py
(somatic precision/recall evaluation -- this ignores the GT and just checks for
presence of alleles). Other tools are qfy.py (which just executes the quantification
step of the analysis pipeline, this requires a
[GA4GH-intermediate](https://github.com/ga4gh/benchmarking-tools/) VCF file), and
[pre.py](doc/normalisation.md), which is hap.py's input cleaning and
variant normalisation step.

Here are some small example command lines. Advanced features like confident call
 / ambiguity / FP regions are also available, see the documentation for each
 tool for these.

Below, we assume that the code has been installed to the directory `${HAPPY}`.

### hap.py

See also [doc/happy.md](doc/happy.md).

```bash
$ ${HAPPY}/bin/hap.py  \
      example/happy/PG_NA12878_chr21.vcf.gz \
      example/happy/NA12878_chr21.vcf.gz \
      -f example/happy/PG_Conf_chr21.bed.gz \
      -o test
$ ls test.*
test.metrics.json  test.summary.csv
```

This example compares an example run of GATK 1.6 on NA12878 agains the Platinum
Genomes reference dataset (***Note: this is a fairly old version of GATK, so
don't rely on these particular numbers for competitive comparisons!***).

The summary CSV file contains all high-level metrics:

| Type          |  TRUTH.TOTAL|  QUERY.TOTAL | METRIC.Recall | METRIC.Precision | METRIC.Frac\_NA | TRUTH.TOTAL.TiTv\_ratio | QUERY.TOTAL.TiTv\_ratio | TRUTH.TOTAL.het\_hom\_ratio | QUERY.TOTAL.het\_hom\_ratio|
|---------------|-------------|--------------|---------------|------------------|-----------------|-------------------------|-------------------------|-----------------------------|----------------------------|
|INDEL          |         9124|         9905 |      0.869406 |         0.978441 |        0.194548 |                     NaN |                     NaN |                    1.463852 |                    1.209105|
|SNP            |        52520|        48078 |      0.894478 |         0.998258 |        0.021070 |                2.081002 |                2.082603 |                    1.595621 |                    1.487599|

These numbers tell us the SNP and indel recall of our query VCF against the
truth dataset. See [doc/happy.md](doc/happy.md) for more documentation and some
advice for their interpretation.

### som.py

Som.py is a simple comparison tool based on bcftools. It does not perform genotype or haplotype matching.

See [doc/sompy.md](doc/sompy.md) for more documentation.

```
${HAPPY}/bin/som.py example/sompy/PG_admix_truth_snvs.vcf.gz \
                    example/sompy/strelka_admix_snvs.vcf.gz \
                    -f example/sompy/FP_admix.bed.gz \
                    -o test
[...]
      type  total.truth  total.query     tp     fp   fn    unk  ambi    recall   recall2  precision        na  ambiguous
1     SNVs        16235        47530  15573  14698  662  17259     0  0.959224  0.959224   0.514453  0.363118          0
3  records        16235        47779  15573  14737  662  17469     0  0.959224  0.959224   0.513791  0.365621          0

ls test.*
test.stats.csv
```

The most relevant metrics are recall and precision. UNK calls are the calls that are outside the
coverage of the truthset.

## Using Python 3 hap.py in variantbenchmarking

The application version is `modern-20261101`; the local image name is
`hap.py-py3`. Override the Nextflow task containers explicitly: changing the host
PATH does not change the container used by the workflow. Run the following
examples from the variantbenchmarking workflow root, not this source directory.

### Singularity (Linux x86_64)

Place `hap.py-py3.sif` on storage accessible to every execution node. Create
`hap.py-py3-singularity.config` with the actual absolute image path:

```groovy
singularity.enabled = true
singularity.autoMounts = true
docker.enabled = false
process {
    withName: 'HAPPY_(HAPPY|PREPY)' {
        container = '/absolute/path/hap.py-py3.sif'
        containerOptions = '--env TMPDIR=.'
        ext.happy_version = 'modern-20261101'
    }
}
```

`TMPDIR=.` uses the writable task directory instead of the read-only image.
Inputs, the image, and the Nextflow work directory must be visible on all nodes.
Follow site-specific bind-mount rules. Nextflow uses `singularity exec`, not the
image runscript.

Keep the samplesheet, reference, truth set, BED, and analysis parameters required
by your analysis. Load the override after any site-specific configurations:

```bash
nextflow run . -profile singularity \
  -params-file params.yaml \
  -c hap.py-py3-singularity.config \
  --analysis germline --variant_type small --method happy \
  --outdir results-hap.py-py3 \
  -work-dir work-hap.py-py3 \
  -with-trace trace-hap.py-py3.txt
```

Supply your own `params.yaml`; it is not provided by the image. Replace existing
output/work directory options instead of duplicating them. Use independent
working directories for initial validation rather than reusing old task results.

### Docker alternative

Import the image archive and add a local alias compatible with the workflow's
registry settings:

```bash
docker load -i hap.py-py3.tar.gz
docker tag hap.py-py3 quay.io/local/hap.py-py3:latest
```

The `quay.io/local` alias is local only; it does not imply publication. Create
`hap.py-py3-docker.config`:

```groovy
docker.enabled = true
singularity.enabled = false
process {
    withName: 'HAPPY_(HAPPY|PREPY)' {
        container = 'quay.io/local/hap.py-py3:latest'
        containerOptions = '--pull=never --env TMPDIR=.'
        ext.happy_version = 'modern-20261101'
    }
}
```

Use the preceding Nextflow example with `-profile docker` and
`-c hap.py-py3-docker.config`. Preload the image on every execution node. Other
workflow tools retain their original images. Do not combine this override with
an older development-image override.

Only `HAPPY_HAPPY` and `HAPPY_PREPY` are replaced, not som.py or standalone RTG
modules. Engine selection still comes from workflow arguments. To request
vcfeval, add the following selector inside the same `process` block, retaining
other arguments required by your analysis:

```groovy
    withName: 'HAPPY_HAPPY' {
        ext.args = '--engine=vcfeval --leftshift'
    }
```

Inspect task `.command.run` files and version records to confirm the actual
image and application version. `HAPPY_PREPY` runs only when the selected
preprocessing mode requires it. Targeted three-engine image acceptance does not
establish full workflow or historical golden-suite equivalence. Validate your
own inputs, references, and stratifications. See the
[container build and usage guide](doc/modern-build.md).

## Installation

The modernization build uses **CPython 3.14**, an explicitly provisioned dependency
prefix, and version-independent private modules in `prefix/lib/hap.py`.
Configuration does not download software or execute Git.

### Docker (recommended)

The rebuilt `hap.py-py3` image targets Linux x86_64. It has not been published
to a public registry. Transfer it with `docker save` / `docker load` before use;
do not attempt to pull this local tag. See the
[container build and usage guide](doc/modern-build.md) for build instructions,
mounts, and validation limitations.

```bash
mkdir -p results/tmp
docker run --rm --pull never --network none \
  --user "$(id -u):$(id -g)" --workdir /out \
  -e TMPDIR=/out/tmp -e RTG_MEM=4G \
  -v "$PWD/data:/data:ro" -v "$PWD/results:/out" \
  hap.py-py3 \
  hap.py /data/truth.vcf.gz /data/query.vcf.gz \
  -r /data/reference.fa -f /data/confident.bed \
  -o /out/benchmark --engine vcfeval --threads 4 --scratch-prefix /out/tmp
```

The `data` directory must contain a matching reference FASTA/FAI, truth and query
BGZF VCFs with indexes, and a confident-region BED. Known historical golden-suite
differences remain; successful builds and targeted acceptance do not imply that
all historical tests pass.

### Helper script

`external/bootstrap.sh` explicitly provisions dependencies; `install.py` installs
to a new versioned directory. Configuration performs no downloads or implicit
dependency installation. Preserve absolute dependency prefixes when installing
on the host: moving the installation directory alone is insufficient. See the
container guide for details.

### Compiling from source with CMake

```bash
cmake -S . -B .work/build   -DHAPPY_DEPENDENCY_PREFIX="$PWD/.pixi/envs/default"   -DPython3_EXECUTABLE=/absolute/path/to/python3.14   -DHAPLOTYPES_VERSION=modernization -DCMAKE_BUILD_TYPE=Release
cmake --build .work/build --parallel 16
cmake --install .work/build --prefix "$PWD/install"
```

For an independent instrumented build use another build directory and
`-DCMAKE_BUILD_TYPE=Debug -DHAPPY_SANITIZERS=ON`. Never reuse a noninstrumented
build cache for sanitizer acceptance.

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
