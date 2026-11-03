# Python runtime structure

The runtime uses functions for orchestration and scientific transformations,
existing polymorphism for caller features and ROC strategies, and small objects
only where state or resource ownership warrants them.

## Runtime modules

- `hap.py`: `build_parser` declares the CLI; `preprocess_inputs` preserves the
  truth-before-query dependency for gender and locations; `compare_inputs`
  dispatches native engines and returns an immutable `ComparisonResult`.
- `pre.py`: standalone and embedded VCF preprocessing. Its public `preprocess`
  interface remains compatible.
- `qfy.py`: `resolve_regions` resolves confident and stratification files in
  insertion order; `quantify` runs quantification and writes reports.
- `som.py`: existing pure summary/extended transformations remain functions;
  `coverage_normalization` computes mean per-chromosome BAM coverage times three.
- `Haplo`: native engine adapters, quantification, and ROC transformations.
  Engine execution models differ; no common engine inheritance tree is imposed.
- `Somatic.FeatureSet` and `Tools.roc.ROC`: existing registries and caller-specific
  strategies remain unchanged. Caller field order and transformations are not
  merged merely to reduce line count. A shared pure list-value formatter removes
  duplicate conversion logic without adding inheritance.
- `Tools.scratch.ScratchWorkspace`: cleans only explicitly registered artifacts.
  Keep mode preserves them. Missing files are harmless; other filesystem errors
  are logged without hiding the original comparison exception. Register only
  paths owned by the current run. Directory removal requires explicit opt-in;
  symlinks are removed without traversing their targets.

## Compatibility and execution

Existing CLI arguments, defaults, public helper names, output schemas, native
algorithms, and top-level multiprocessing workers are retained. Comparison
results carry paths, format, and ROC mapping, not pools or open file handles.
Preprocessing still mutates the legacy argument namespace; this is not a fully
immutable pipeline. Pools close and join after successful reporting; failed runs terminate and join.

CMake copies the complete Tools/Haplo/Somatic directories, including the new
resource module. Installed-layout integration tests exercise these copied files,
not the older Python source embedded in the container image.

## Verification

Use the complete dependency environment for unittest discovery. Also run
`test_native_engines.py NATIVE_BIN RTG_EXECUTABLE` separately; it is not executed
by unittest discovery. Cover engines with one/two threads, ROC/UNK, spawn and
forkserver, and run `test_scmp_distance_blocks.py SCMP_EXECUTABLE` separately.
Compare pre/post-refactor scientific tables without regenerating expected data.

Small-fixture timing is a regression signal, not a whole-genome performance
claim. Report repeated runs and variability. Linux wait4 maximum RSS is not the
sum of all concurrently running process memory.

## Remaining work

This is an incremental refactor, not full runtime compatibility update. Legacy broad
exception handling, long caller extraction functions, mutable argparse state
remain outside this targeted refactor. Quantification now uses literal argv; RTG commands execute without a shell. No wholesale
formatting, dependency upgrade, engine hierarchy, or golden-data refresh is
included. The vcfeval preserve-info path now annotates sample-free source INFO onto
comparison records and indexes the result. Native tests assert source INFO,
comparison sample names, counts, ROC and one/two-thread equivalence. Historical full-reference suites have not been established by the
small-fixture tests.

## Invalid-record quarantine

`Tools.invalid.split_invalid` owns streaming classification and durable quarantine
publication behind one interface. `pre.preprocess` owns its normal-record
scratch file and removes it on success or failure. Both CLI entry points use this
same seam before `vcfcheck`; the comparison engines retain their strict native
GT/AD safety checks. Raw GT text is inspected because pysam's decoded GT tuple
maps an out-of-range allele index to `None`, losing the distinction from a true
missing genotype. BCF records are serialized via HTSlib before the same checks.
