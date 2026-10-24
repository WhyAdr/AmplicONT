# FAPROTAX2 production-path validation

**Date:** 2026-10-24
**Validated mode:** taxonomy
**Sequence mode:** blocked; not implemented or validated

This note records the native-Windows integration evidence for FAPROTAX2. The successful pipeline input is the repository's public synthetic `synthetic_minimap2` fixture. It verifies the adapter, runner, manifests, outputs, and preservation of the legacy module; it is not a result from the user's new 16S datasets.

## Full pipeline run

The completed run is recorded in `output/faprotax2-taxonomy-example/run_manifest.json`. The manifest binds it to candidate commit `9e561cc39dd072c244b04bee56ff9f9141896c4c`, records `git_dirty: false`, a synchronized lock, schema revision 6, and SHA-256 entries for 49 source files. It owns 88 files overall, including 13 FAPROTAX2 outputs.

The single-sample run requested all nine modules. Six completed: QC, alpha diversity, composition, kreport, legacy FAPROTAX, and FAPROTAX2. Beta diversity, ordination, and shared taxa were recorded as structured skips because the fixture has one sample. The pipeline completed transactionally and passed both the independent FAPROTAX2 run verifier and the release integration verifier. The separate `--validate-only` run passed and reported zero filesystem mutations.

The FAPROTAX2 input had one positive-count eligible prokaryotic feature. Its mapper output contained 101 functions, with 7 positive and 94 absent feature-function states and no NA states in this fixture. Coverage reported 5 total input alignments, 1 eligible taxon alignment, and 4 excluded input alignments. This one-feature result is a contract check only; the nine-module run does not test multi-sample comparative analyses.

The existing FAPROTAX 1.2.12 module also completed in the same invocation. A separate legacy baseline run completed and passed release verification; the legacy module's test file passed 4 tests and 27 expectations. The FAPROTAX2 work leaves the legacy module and its output directory intact.

A second single-sample pipeline run used the maintained zero-classified synthetic fixture. FAPROTAX2 completed as a structured skip with `No positive-count Bacteria or Archaea rows`; the independent run helper verified that skip, and the fixture inputs remained unchanged. This confirms the eligibility skip does not require launching the mapper for that input.

## FAPROTAX2 contract and interpretation

The adapter ran the audited, unmodified FAPROTAX-mapper 2.0 beta script. Its provenance records the official archive SHA-256 and exact hashes for the mapper and six reference assets. These upstream files and mapper-specific libraries are local runtime assets; the repository does not bundle them or download them automatically.

The module emits these 13 files under `09_FAPROTAX2/`:

- `faprotax2_function_abundance.tsv`, `faprotax2_coverage.tsv`, and `faprotax2_feature_function_assignments.tsv`
- `faprotax2_top_functions.png` and `faprotax2_provenance.json`
- The mapper input, identity crosswalk, raw function/state/probability matrices, report, stdout, and stderr

Function abundance is calculated from binary predicted-presence states and the input taxon counts. NA states contribute zero to that calculation; probabilities are retained as unitless annotations and do not weight abundance. Function assignments overlap, so their summed abundance can exceed the input total and is not normalized to 100%. The coverage table separates total input, excluded input, eligible input, predicted presence, known absence, and unknown-without-presence. Taxonomy translation omits the NCBI kingdom rank, retains the source taxon path, and records the seven mapper ranks.

The synthetic pipeline run had integer counts, as required by the existing pipeline input contract. Separate synthetic matrix tests confirm that fractional arithmetic is preserved rather than rounded; they do not add a fractional-count input contract or verify the mapper's probability-weighted mode. Invalid state/probability combinations are rejected, matrices must match by feature/function/sample IDs and order, and missing or altered pinned assets fail validation. NA is accepted by the adapter's matrix contract, but no NA states occurred in the known-lineage native mapper fixtures, so a positive NA-producing case has not been empirically observed.

An earlier native-Windows mapper POC ran the three known-lineage fixture twice without a random seed. The four raw result tables were byte-identical between those two runs. This is a local repeatability observation, not a cross-platform determinism guarantee. More details and hashes are in [`faprotax2-poc-results.md`](faprotax2-poc-results.md).

The pipeline seed is not forwarded to the external mapper; no upstream seed flag
was invented. Reproducibility claims are limited to the observed pinned runs.

## Native runtime and timeout evidence

The production-path run used Windows 11, R 4.5.3, the synchronized project `renv` library, one mapper thread, and `timeout_seconds: 600`. Provenance recorded 27 resolved child package version and DESCRIPTION hashes. The optional mapper library supplied `castor` 1.8.7 and `naturalsort` 0.1.3; other dependencies resolved from the project library and base R. This is a recorded library overlay, not a hard-isolated R installation. No package installation occurred as part of the mapper run.

The mapper timeout is measured in seconds. An earlier implementation multiplied this value by 1,000 based on the wrong unit assumption; the mistake was corrected after checking `processx::run`. A real child-process regression now runs a fake mapper that sleeps for 30 seconds with a one-second timeout. It verifies an explicit timeout error, child cleanup before the completion marker, and no function summary publication. The actual 1-second pipeline timeout trial also recorded the timeout diagnostic, aborted the transaction, and preserved all 16 previously published files byte-for-byte. Its log additionally shows the staging-output census diagnostic for intermediate files left in the failed stage; those files were discarded and not published.

## Test and verification results

The final native Windows R test suite completed with 232 tests, 1,494 passing expectations, 14 warnings, zero failures, and zero skips (909.4 seconds). The warnings came from existing cohort tests, including vegan test-data conditions and UpSetR/ggplot2 deprecation notices; the suite passed.

Additional checks:

- Focused FAPROTAX2 tests: 112 expectations passed, including the real `processx` timeout regression.
- Native Windows Python suite: 85 tests ran in 23.610 seconds: 83 passed and two symlink tests were skipped because the account lacked Windows symlink privilege. A separate sandboxed attempt was inconclusive and is not counted.
- Local Node 18.8.0 test: one TAP file passed its client assertions; CI separately uses Node 24.11.1.
- Release integration verifier: passed for `faprotax2_taxonomy_contract`.
- Legacy FAPROTAX regression tests: 4 tests and 27 expectations passed.

The final suite output and run evidence are retained locally under `output/faprotax2-audit/runtime/`; those logs and mapper assets are ignored and are not distributed with the repository.

## Reproducing taxonomy mode

The public fixture config is [`tests/fixtures/faprotax2/config.taxonomy.yml`](../tests/fixtures/faprotax2/config.taxonomy.yml). It uses synthetic inputs, enables taxonomy mode, and points to the local audited mapper assets and mapper dependency library under `output/`. A clean checkout must first provision those assets and packages at the configured paths, or edit the explicit paths to their local locations. The mapper does not fetch them automatically.

From the repository root in native PowerShell, choose a fresh output directory and run:

```powershell
$env:RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE'
$lockedLibrary = (Resolve-Path 'renv/library/windows/R-4.5/x86_64-w64-mingw32').Path
$env:R_LIBS = $lockedLibrary
$env:LC_ALL = 'English_United States.utf8'

Rscript --vanilla analysis/00_run_pipeline.R `
  --config tests/fixtures/faprotax2/config.taxonomy.yml `
  --modules faprotax2 `
  --output-dir 'output/faprotax2-taxonomy-only-repro'
```

Use `--overwrite` only if the selected output directory is one you intend to replace. The tested nine-module invocation used the same config with `--modules qc,alpha,beta,composition,ordination,shared,kreport,faprotax,faprotax2` and produced the run described above.

## Sequence-mode limit

Setting `faprotax2.mode: sequence` is explicitly rejected by the current runtime gate. The one-feature native trial reached APPLES placement output but the mapper rejected the resulting tree because its topology differed from the FAPROTAX reference tree; no functional abundance output was produced. Sequence mode therefore has no verified usage command. See [`faprotax2-native-sequence-trial.md`](faprotax2-native-sequence-trial.md) for the complete evidence. Do not interpret the partial placement as a successful FAPROTAX2 sequence run.
