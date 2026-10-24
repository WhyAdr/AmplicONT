# FAPROTAX2 migration assessment

Phase 1 was investigated on 2026-10-24 in native Windows PowerShell. The live
checkout is AmplicONT, not the former repository name printed in the submitted
plan. Starting `HEAD` and fetched `origin/main` both equal
`42572d15b059ef29efe46ef5fe67010d06b1756a` (VERSION 0.5.0). The only initial
untracked work was the submitted migration plan and private input directory;
both are preserved. The submitted plan remains a local, untracked source.

## Native runtime and unchanged legacy baseline

| Component | Verified result |
|---|---|
| OS / shell | Windows 11 x64 build 26200 / PowerShell 7.6.5 |
| R | 4.5.3 ucrt, native Windows executable |
| Python | 3.12.10, native Windows executable |
| Existing library | Repository `renv/library/windows/R-4.5/x86_64-w64-mingw32` |
| Lock comparison | All 69 locked packages installed at exactly the locked versions |
| Legacy inference | microeco 2.3.0; embedded FAPROTAX 1.2.12 |
| Process adapter | processx 3.9.0 |
| Optional mapper packages | optparse 1.8.2 and stringr 1.6.0 present; castor missing |
| CRAN Windows binary query | castor 1.8.7 and naturalsort 0.1.3 available for R 4.5 |
| Laptop resources | 8 logical CPUs; 8.28 GB total RAM; about 1.23 GB free at inspection; about 172 GB disk free |

The initial restricted-process probe incorrectly reported that installed packages
were absent. Running the same native R probe outside that restriction recognized
the installed library. This was an execution sandbox issue, not a broken renv
installation. The working library was not restored, upgraded, or overwritten.

Environment for all baseline R commands:

```powershell
$env:RENV_CONFIG_AUTOLOADER_ENABLED = 'FALSE'
$env:R_LIBS = (Resolve-Path 'renv/library/windows/R-4.5/x86_64-w64-mingw32').Path
$env:LC_ALL = 'English_United States.utf8'
```

The actual prior Bioslurry run record dated 2026-10-22 was read locally and its
environment and eight-module invocation were checked against the current CLI.
The new baseline used the repository's existing config and inputs, with a fresh
output root. No biological inputs were transferred to any service.

```powershell
Rscript --vanilla -e 'r <- testthat::test_dir("tests/testthat", filter="^faprotax$", reporter="summary"); d <- as.data.frame(r); if(any(d$failed > 0 | d$error)) quit(status=1)'
Rscript --vanilla analysis/00_run_pipeline.R --config config.yml --modules qc,alpha,beta,composition,ordination,shared,kreport,faprotax --output-dir 'output/faprotax2-audit/legacy baseline Rényi'
Rscript --vanilla tests/verify_release_run.R 'output/faprotax2-audit/legacy baseline Rényi'
```

All three commands exited zero. Legacy tests: **4 tests, 27 expectations, zero
failures, errors, or skips**. The pipeline completed transactionally and its
independent release verifier passed. QC, alpha, composition, kreport and legacy
FAPROTAX completed; beta, ordination and shared executed as structured
single-sample skips. The manifest inventories 48 maintained source files and 71
owned outputs. The taxonomy diagnostics retained 46 unresolved nodes and 26
conflict records; successful execution does not resolve those biological labels.
The output path contains both spaces and Unicode; the UTF-8 locale also preserved
the `Rényi` probe. Logs, package comparison, manifest and artifacts remain local
under ignored `output/faprotax2-audit/`.

## Upstream and inputs

See [official contract and archive inventory](faprotax2-upstream-contract.md),
[source integration audit](faprotax2-architecture-assessment.md), and
[private-input assessment](faprotax2-input-assessment.md).

The downloaded **mapper**, rather than FAPROTAX-db, is 12,788,883 bytes with
SHA-256 `859b99133b55f1e4f637e6ac14a804f3321ed27f76ee03d95e3f215382206bac`.
Its embedded BSD-style conditions permit use/redistribution subject to retained
notices and disclaimers and no author endorsement. No upstream code, reference
archive or derived biological reference fixture is bundled in this phase.

Six real 16S FASTA/count-table pairs have exact feature correspondence. Their
metadata identifies Illumina/DADA2, so an end-to-end run on them cannot be called
an ONT validation. Original files remain unchanged and ignored. Private absolute
paths, sample/feature labels, sequence bases and inventory hashes are kept out of
the published reports.

## Gate P1 and next gate

The unchanged native legacy baseline passes. Official source and native `--help`
establish separate taxonomy and representative-16S input interfaces. Taxonomy is
supported by the official contract; native prediction must still pass Phase 2.
Sequence inference is supported upstream but needs MAFFT plus EPA-NG or APPLES
and contains shell-based subprocess operations; native end-to-end feasibility
is not established. Neither mode is implemented or empirically validated by
this documentation commit. No syntax or prediction result has been inferred
from a third-party wrapper.

Phase 2 must install optional mapper dependencies into an isolated native library,
preserve the existing lock, execute each mode independently, examine raw numeric
outputs, and test offline reruns and controlled failures. Production integration
is permitted only for modes that pass that gate. Conservative trials use one
thread because available laptop memory is limited.
