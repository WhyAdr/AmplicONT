# FAPROTAX backend validation: taxonomy comparison

**Date:** 2026-10-24

**Status:** Two-backend taxonomy comparison completed; three-way comparison blocked

**Run revision:** `1778e3a6370c019165d2dc6a65356610f0b451a7` (manifest records `git_dirty: false` for maintained source files; schema revision 6)

This report records an independent comparison of the existing FAPROTAX 1.2.12 / microeco route and FAPROTAX-mapper 2.0 beta in taxonomy mode. Both modules completed in the same native Windows run with the same `minimap2_alignments` input count basis. The run used one sample. No sequence-mode result is included.

## Independent checks

The comparison harness reconstructs eligible Bacteria and Archaea rows from the original abundance table, requires the same sample order and `Taxon_######` IDs in the FAPROTAX2 input, verifies the source-taxonomy crosswalk, and compares FAPROTAX2 matrix products to the mapper summary. For legacy FAPROTAX, it reconstructs the feature-by-function binary matrix from the published assignments and multiplies it by the original eligible count matrix. It then checks the legacy function table, percentages, and coverage against those independent values.

All **1,836** eligible source taxa and their counts matched exactly across the original input and the FAPROTAX2 mapper input. Both modules reported `minimap2_alignments`. FAPROTAX2 recorded `normalize: none` and `probability_weighting: false`; its probabilities were retained as annotations and did not weight the function abundance. The comparison rechecked all **12** file fingerprints recorded for input and FAPROTAX2 resources. The pinned official archive and mapper source hashes matched the audited release.

The raw FAPROTAX2 state matrix contained 185,436 cells: **15,524 present, 169,614 absent, and 298 unknown (NA)**. All 1,836 eligible taxa had at least one positive function state; there were no all-zero taxa or taxa with unknown states and no positive function state. That is output-level predicted-function coverage. It does not establish successful or accurate placement, and the adapter warns that unrecognized taxonomies may receive broad fallback predictions.

Legacy FAPROTAX had at least one positive assignment for **924 of 1,836** eligible taxa, corresponding to **52,403** eligible counts. FAPROTAX2 had a positive function state for all 1,836 taxa, corresponding to **80,556** eligible counts. These backend outputs differ in their assignment results; the comparison does not treat either count as an accuracy estimate.

## Observed function labels and abundance

The FAPROTAX2 configured state matrix contained **101** function labels, of which **72** had positive abundance. The legacy run emitted **78** labels with observed assignments. Exact label names overlapped for **77** labels; one legacy observed label, `denitrification`, was absent from the FAPROTAX2 configured catalog. Twenty-four FAPROTAX2 catalog labels were absent from the legacy run's observed assignments. This describes the two per-run outputs; the legacy output does not enumerate definitions with no assignments, so this is not a comparison of the complete embedded database catalogs.

Of the observed positive labels, **61** were positive in both outputs, **17** were positive only in the legacy output, and **11** were positive only in FAPROTAX2. The legacy-only labels were:

`aerobic_anoxygenic_phototrophy`, `aliphatic_non_methane_hydrocarbon_degradation`, `anoxygenic_photoautotrophy_H2_oxidizing`, `anoxygenic_photoautotrophy_S_oxidizing`, `chitinolysis`, `chlorate_reducers`, `dark_oxidation_of_sulfur_compounds`, `dark_sulfite_oxidation`, `denitrification`, `dissimilatory_arsenite_oxidation`, `human_pathogens_all`, `human_pathogens_diarrhea`, `human_pathogens_meningitis`, `human_pathogens_nosocomia`, `human_pathogens_septicemia`, `invertebrate_parasites`, and `plastic_degradation`.

The FAPROTAX2-only positive labels were:

`acetate_oxidation`, `alcohol_fermentation`, `collagen_degradation`, `facultative_anaerobes`, `gelatin_degradation`, `iron_oxidation`, `lactate_oxidation`, `lactic_acid_fermentation`, `obligate_anaerobes`, `pectin_degradation`, and `sulfur_reduction`.

For the 77 exact-name labels present in both backends' run outputs, **5** sample/function counts matched and **72** differed. Identical labels do not establish equivalent definitions, taxonomic rules, or validated biological meaning. Labels such as `human_pathogens_*` are database phenotype associations from taxonomic assignment; they do not establish pathogen genes, strain-level identity, pathogenicity, virulence, clinical risk, expression, or disease. Detailed sample-, feature-, and function-level tables are kept only in ignored local `output/` evidence; the report contains no sample IDs, feature IDs, raw taxonomy rows, DNA, private paths, or input hashes.

Functions overlap within each backend. The sum of legacy function counts was **158,498**, and the sum of FAPROTAX2 function abundances was **694,627**, both greater than the **80,556** eligible input count. These sums are not partitions; do not force overlapping function categories into mutually exclusive 100%. FAPROTAX2's predicted-presence abundance of 80,556 counts each eligible input taxon once if it has any positive function state; that coverage value is not evidence of placement accuracy.

## Software provenance

- Legacy: microeco **2.3.0**, embedding FAPROTAX **1.2.12**. The installed microeco `DESCRIPTION` file SHA-256 was `8c6aab860d75b5383dee5393ce9d624706add53259d7a8aaa52dca8005022b92`. This is a package metadata-file fingerprint, not a hash of the complete installed package tree.
- FAPROTAX2 official mapper archive SHA-256: `859b99133b55f1e4f637e6ac14a804f3321ed27f76ee03d95e3f215382206bac`.
- FAPROTAX2 official mapper script SHA-256: `8bd6ff97aa53f17fb0122b98c436275c230e34629f7a83c446b2b3854e53a64f`.

The run manifest was bound to the clean commit above. The mapper and reference data were local copies at paths containing spaces and Unicode; their release hashes were unchanged. No upstream assets are bundled by this report.

## Validation evidence and limits

- The native full R suite passed **232 tests / 1,494 expectations**, with **14 warnings**, **0 failures**, and **0 skips**. Warnings were from existing cohort-test conditions and dependency deprecations.
- Focused FAPROTAX2 tests passed **112 expectations**, including the real `processx` timeout regression. Legacy FAPROTAX regression tests passed **4 tests / 27 expectations**.
- The native Windows Python suite ran **85 tests**: **83 passed** and **2 symlink tests skipped** because the account lacked Windows symlink privilege. The local Node 18.8.0 client test passed; release CI uses Node 24.11.1.
- The independent FAPROTAX2 run validator passed on the completed real-data run. The P5 comparator independently validated both function-abundance products and coverage arithmetic, sample order, feature IDs, and original input/resource fingerprints.
- The comparator rejected a private evidence target under the completed run root before creating it; a positive rerun with the evidence target outside that root passed.
- The pinned unmodified FAPROTAX2 taxonomy mapper produced byte-identical raw tables in two earlier unseeded synthetic repeats. This is a repeatability observation for those fixtures, not a general determinism guarantee.
- An actual one-second mapper timeout aborted transactionally and preserved all 16 previously published files byte-for-byte. The no-prokaryote path also passed as a structured skip without mapper computation.

Sequence mode remains blocked from P2: the native trial generated an APPLES jplace result but FAPROTAX rejected its re-estimated tree for reference-tree topology mismatch before annotation. P4 sequence implementation is not complete or validated. Therefore this report is **not** a three-way comparison, makes no sequence-mode success claim, and does not compare placed sequences, rank concordance, or biological accuracy. See [the production validation record](faprotax2-production-validation.md) and [the native sequence trial](faprotax2-native-sequence-trial.md).

## Reproduce the comparison

The harness reads an already completed run and writes all detailed evidence beneath the ignored repository `output/` directory. It does not run either mapper. From a native Windows checkout with the locked R packages available:

```powershell
$env:RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE'
$env:R_LIBS = (Resolve-Path 'renv/library/windows/R-4.5/x86_64-w64-mingw32').Path
$env:LC_ALL = 'English_United States.utf8'
Rscript --vanilla tests/compare_faprotax_backends.R `
  '<completed-run-directory>' `
  'output/private-faprotax-comparison-evidence'
```

The first path must contain `run_manifest.json` and completed outputs from both modules. The second path must remain under `output/`; the harness refuses to write private per-sample and per-feature evidence elsewhere.
