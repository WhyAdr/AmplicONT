# FAPROTAX-mapper 2.0 taxonomy POC results

**Date:** 2026-10-24
**Scope:** Native Windows taxonomy route only; official mapper source was not modified.

## Runtime and inputs

- Windows R 4.5.3 (2026-03-11 ucrt); Python 3.12.10 was available but not used by the mapper POC.
- Isolated mapper library: castor 1.8.7 and naturalsort 0.1.3. Locked library supplied optparse 1.8.2 and stringr 1.6.0. Child R reported the isolated mapper library first in its library paths. No package installation occurred during the POC.
- Official archive: FAPROTAX-mapper_v2.0.tar.gz, 12,788,883 bytes, SHA-256 859b99133b55f1e4f637e6ac14a804f3321ed27f76ee03d95e3f215382206bac.
- Official script SHA-256: 8bd6ff97aa53f17fb0122b98c436275c230e34629f7a83c446b2b3854e53a64f.
- Synthetic input SHA-256: 491badafded5070c9916c89a2bdcfc19e7684fd9891cacbfbbf91e85bb24dd63.
- Fixture: taxonomy_known_species.tsv, three reference-matching lineages (Nitrosomonas europaea, Pseudomonas aeruginosa PA38182, Bacillus subtilis), two samples (Rényi_α, Control_β), integer totals 53 and 51. Fixture IDs are unique. No mapper/reference data or user sequences were added to the repository.

## Successful taxonomy runs

The unmodified mapper was launched as a native Rscript.exe process with processx argument vectors, a child R profile that places the isolated mapper library first, explicit --Nthreads 1, --normalize none, --hsp_algorithm config, and --unknown_ref_states config. The default per-function configuration yielded 101 function rows. Output directories and filenames contained both spaces and Unicode; both runs exited 0.

| Run | Elapsed | Sampled peak RSS | State / probability / function / overlap table SHA-256 |
|---|---:|---:|---|
| Run one | 70.764 s | 352,776,192 B | all table hashes below |
| Rényi repeat two | 71.238 s | 353,120,256 B | identical to run one |

The observed state matrix had 22 present, 281 absent, and 0 NA cells across the 3 × 101 feature-function pairs. All 303 annotation probabilities were finite and in [0, 1] (observed min 0, max 1); every observed 0/1 state agreed with the configured thresholds. All feature IDs and both sample labels were retained. The harness independently recalculated the function table as the transpose of the 0/1/NA state matrix, with NA treated as zero, multiplied by the input feature-by-sample count matrix. The numeric matrices matched exactly at tolerance 0.

Per-feature presence counts were 7 for feature_neuropaea, 11 for feature_aeruginosa, and 4 for feature_subtilis. Summed function abundances were 379 and 389, while input totals were 53 and 51. This is expected overlap among functions, not a partition of sample counts. Probability weighting was not enabled.

The two unseeded runs produced byte-identical matrix table files even though the scratch paths changed:

- Function abundance: 7a1352113ed41fa4cdb4faaa9573f055bf651ca0df7c4f670ce3adb61583f661
- Feature-function states: 562de3e569199aea77a308f2b5024192f8ad8bee6d98d3996b14cfadbd824731
- Annotation probabilities: 6697a74d369327d91d9cfb29cd9bf4ec57200f74c854ea21dead11d72b45b92c
- Function overlaps: b1eef0dacb522cdbbd466ae91b2d05f9eb321fb04d9fc4064d18a662af772804

The report and mapper log are retained separately for each run under ignored output scratch. Their timestamps and recorded command paths naturally differ. The repeat used local mapper/reference files and no network inputs or installs; operating-system network access was not disabled.

## Boundary and failure observations

- The synthetic unknown lineage in taxonomy_unknown_root.tsv completed with exit 0 under --placement_failure error, and the report retained one processed feature. The reproducible harness selects one configured function for this boundary probe and does not assert its predicted state. An earlier exploratory unfiltered run on the same taxonomy and counts (before the fixture ID was renamed) completed all 101 functions and reported 5 feature-function associations. Unknown taxonomy text therefore cannot be equated with a placement failure or with zero function assignments.
- taxonomy_invalid_count.tsv returned exit 1 before writing a function table. The official diagnostic states that a sample column was nonnumeric and recommends --non_numeric_columns omit or sample filters.
- The harness killed and reaped a separate 30-second sleeping child Rscript before its completion marker was written. This verifies child cleanup logic; it is not a forced-interruption trial of a running FAPROTAX mapper.

## Reproduce

Run from the repository root in native PowerShell. Use a fresh output directory whose path includes a space and Rényi.

~~~powershell
$env:RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE'
$mapperLibrary = (Resolve-Path 'output/faprotax2-runtime/library').Path
$lockedLibrary = (Resolve-Path 'renv/library/windows/R-4.5/x86_64-w64-mingw32').Path
$env:R_LIBS = "$mapperLibrary;$lockedLibrary"
$env:LC_ALL = 'English_United States.utf8'
Rscript --no-save --no-restore tests/run_faprotax2_poc.R --mapper 'output/faprotax2-audit/upstream/unpacked/FAPROTAX-mapper_v2.0/faprotax.R' --reference-dir 'output/faprotax2-audit/upstream/unpacked/FAPROTAX-mapper_v2.0/FAPROTAX_ref_data' --mapper-library $mapperLibrary --locked-library $lockedLibrary --output-dir 'output/faprotax2-audit/new taxonomy POC Rényi'
~~~

The harness saves arguments, child stdout/stderr, mapper log/report, resource samples, hashes, and summary under the ignored output directory. It uses the official table contract with -t column:taxonomy; it does not implement a production adapter.

An additional lead-agent rerun completed all 101 functions with R's `download.file`,
`install.packages`, and `url` entry points traced to raise an error. Its function,
state and probability tables matched the corresponding hashes above. This verifies
cached operation with those R entry points blocked; it does not claim an OS
firewall or complete network isolation. The tested guard is preserved in
`tests/faprotax2_offline_bootstrap.R`; set `FAPROTAX2_MAPPER` to the official script
and pass the same mapper arguments after that bootstrap script. Local evidence is
under `output/faprotax2-audit/offline taxonomy Rényi/`.

## Limits

- The successful known-taxon run produced no NA state cells. The harness permits NA and verifies the documented NA-to-zero multiplication rule if such cells occur, but an actual upstream NA state was not observed in this fixture. A separate threshold-boundary fixture would be needed to claim that case empirically verified.
- The unknown-lineage probe shows that --placement_failure error is not a sufficient taxonomy-coverage check. Any production coverage measure must distinguish accepted fallback placement, function-state assignments, and features with no positive function.
- This phase did not validate sequence alignment/placement, probability-weighted output, or forced termination of the active mapper process.
