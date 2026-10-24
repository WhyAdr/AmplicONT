# FAPROTAX2 migration delivery and reproducibility

This migration delivers the demonstrated native-Windows taxonomy route. It
does not complete the blocked sequence route or a three-backend comparison.
`VERSION` remains 0.5.0, and no release tag is created. The legacy
`analysis/08_faprotax.R`, its output contract, and `renv.lock` are unchanged.

## Phase disposition

| Phase | Result | Evidence |
|---|---|---|
| P1 investigation | Complete | [Audit](faprotax2-migration-assessment.md), [official contract/license](faprotax2-upstream-contract.md), [private-input inventory summary](faprotax2-input-assessment.md) |
| P2 native proof | Taxonomy passed; sequence failed independently | [Mapper proof](faprotax2-poc-results.md), [sequence trial](faprotax2-native-sequence-trial.md) |
| P3 production integration | Taxonomy passed | [Clean nine-module run, numerical oracle, rollback and regression checks](faprotax2-production-validation.md) |
| P4 sequence integration | Blocked; production mode rejected | [Feasibility gate](faprotax2-sequence-feasibility.md) |
| P5 comparison | Two taxonomy backends passed; three-way blocked | [Independent comparison](faprotax2-validation-report.md) |
| P6 documentation/setup | Native setup and offline contracts verified; final pushed CI reported in the delivery handoff | [README tutorial](../README.md#native-windows-setup-and-example) |

## Verified publication history before the final documentation phase

The live starting local and remote `main` was
`42572d15b059ef29efe46ef5fe67010d06b1756a`, later than the plan's baseline.
Every publication used a normal fast-forward push followed by an exact
`git ls-remote origin refs/heads/main` comparison. These are recorded remote
heads at each completed publication, not claims that every intermediate CI run
passed. The final phase's commit/head and CI are supplied in the handoff.

| Change | Commit | Verified pushed remote head |
|---|---|---|
| P1 audit | `f56457f7c70bfebffceffe6a90fb36135c5562c1` | Same SHA |
| Baseline Linux missing-value test correction | `af6b1692fe83649faded28f834a0f98433bb009f` | Same SHA |
| P2 native taxonomy proof | `f610be557117e39050ba43ca4e6831f9e877435c` | Same SHA |
| P3 adapter and compatible readers | `4a17dafcb241ce26e7b4dfadd2b295de885ef604` | Same SHA |
| P3 production registration | `9e561cc39dd072c244b04bee56ff9f9141896c4c` | Published together with the next commit at `ec7818b8cda9f24d62d858d6958afb96feb737cd` |
| P3 independent release gate | `ec7818b8cda9f24d62d858d6958afb96feb737cd` | Same SHA |
| P4 blocked native sequence evidence | `1778e3a6370c019165d2dc6a65356610f0b451a7` | Same SHA |
| P5 paired taxonomy comparison | `71bac2c2eb16ee048593a88a3f5fb0b64c94d939` | Same SHA |
| Windows short-path test correction | `8072bc158b5dfd168af01851ee36f00aa0a91989` | Same SHA |

No PR, force push, release tag, private input upload or upstream asset bundle
was created. The supplied implementation-plan file remains a separately
preserved untracked user file.

## Native setup verification

`tests/setup_faprotax2_windows.R` checks native Windows R 4.5.3 and all 69 raw
package versions from `renv.lock`. In the inspected installation, 63 packages
reside in the project library and six resolve from R's library. It does not
restore or update that baseline. Mapper assets stay under ignored `output/`.
The helper hash-checks the official mapper archive and seven production
resource identities, then installs only the required optional Windows binaries
into an explicit separate library. Existing wrong-version optional packages or
wrong-hash cached archives cause an error rather than silent replacement.

The additional binary pins are:

| Windows R 4.5 binary | SHA-256 |
|---|---|
| `naturalsort_0.1.3.zip` | `aa835431b76b3e979ee2e99ae36e538dce53a3144cac0ecfc0360d1c6fa10259` |
| `castor_1.8.7.zip` | `dbcc5d83d9bcf63cd0b4de9b94d6079553435a42616e01fffdb4d4085c5df847` |
| `RSpectra_0.16-2.zip` | `d3e22c9a67d04d834a772856fa5372f2b1c04a1f8ff6cea32f6e1669bfb211f1` |

Their pinned URLs use
`https://cloud.r-project.org/bin/windows/contrib/4.5/` followed by the listed
filename. CRAN may replace versioned binaries; a future unavailable download
is a setup failure, not permission to silently choose another version. Keep
the hash-matched cache for offline reuse. The helper prints versions and hashes,
never the full process environment.

Fresh extraction was also exercised in a new path containing spaces and
`Rényi`. The official archive has a root-level AppleDouble metadata sidecar;
the helper allows that exact path and rejects unexpected top-level paths.
Extraction occurs in a unique owned temporary sibling. All seven production
resource hashes must pass before promotion to the final cache directory.
A controlled tar failure removed the partial stage and preserved an existing
sibling; a subsequent native tar retry succeeded. Repeating extraction against
the completed tree validated it without rewriting it. A full setup rerun
installed no packages, and an attempted mapper-library target equal to the
locked project library was rejected with its before/after inventory unchanged.

A fresh optional library in a spaced path was populated from those cached
archives. The baseline retained all 69 exact versions, and production dependency
preflight passed. This library includes the pinned RSpectra binary rather than
depending on the previously installed same-version global binary. They are
distinct binary installations and were checked separately.

The official three-feature/two-sample taxonomy proof then passed twice with
that fresh library, spaces and Unicode in paths, all 101 default functions,
22 presence states, 281 absence states and no NA in this small fixture. The
runs took 71.787 and 71.326 seconds; sampled peak RSS was 335,872,000 and
342,286,336 bytes. All four raw table hashes matched between repeats and also
matched the earlier pinned proof. The unknown-root, invalid-numeric-input and
child-termination probes passed. These repeats reused local assets without
installation or network inputs; network access was not disabled at OS level.

Local logs are retained under ignored `output/faprotax2-runtime/` and
`output/faprotax2-audit/`. Synthetic input hashes may be published in the proof;
private sequences, feature tables, inventory paths and raw results remain local.

## Final clean-source usage check

The README's default setup, validate-only, taxonomy run and independent verifier
commands passed at candidate commit
`0a422aaafe7c6da88455276288f13867335e1e47`. The setup rerun installed no packages.
Validate-only reported zero filesystem mutations. The completed taxonomy-only
run took 70.671 seconds, owned 16 files including the 13 FAPROTAX2 artifacts,
and recorded `git_dirty: false`, `locked: true` and manifest revision 6.
The independent verifier passed the one-feature, 101-function, one-sample
output. No development bypass was used.

A separate legacy-only run also completed from that same clean candidate using
the same public fixture config, in 0.8705 seconds with eight owned files
including the unchanged five legacy artifacts. Its manifest also records a
clean source tree, synchronized lock and revision 6. The new output roots are
`output/faprotax2-taxonomy-repro` and `output/faprotax-legacy-repro`.

After the README's environment/setup block, the independently exercised modes
can be invoked from native PowerShell as follows (choose new output roots when
rerunning):

```powershell
Rscript --vanilla analysis/00_run_pipeline.R `
  --config tests/fixtures/faprotax2/config.taxonomy.yml `
  --modules faprotax2 --output-dir 'output/faprotax2-taxonomy-repro'
Rscript --vanilla tests/verify_faprotax2_run.R 'output/faprotax2-taxonomy-repro'

Rscript --vanilla analysis/00_run_pipeline.R `
  --config tests/fixtures/faprotax2/config.taxonomy.yml `
  --modules faprotax --output-dir 'output/faprotax-legacy-repro'
```

## Regression and CI scope

The complete native R suite passed 232 tests and 1,494 expectations, with 14
warnings, zero failures and zero skips. The focused final PowerShell entrypoint
passed 14 FAPROTAX2 tests and 112 expectations. The native Python suite passed
83 tests and skipped two Windows symlink-privilege cases; the local Node client
test passed. Source parsing, dependency checks, release metadata and whitespace
checks passed. Exact logs and run conditions are in the production report.

The pushed Windows jobs exposed a new test-only path comparison defect:
nonexistent asset suffixes were appended to a short-form Windows temp path,
while the config loader had already resolved the existing ancestor's long name.
The test now canonicalizes that existing ancestor before comparing exact
suffixes. The three failed assertions did not indicate an incorrect mapper
path. The correction does not weaken input, resource or matrix validation.

CI restores the existing lock on Ubuntu and Windows. It runs the adapter
contracts offline and adds an explicit native PowerShell invocation on Windows;
the full licensed mapper is validated locally rather than downloaded in CI.
The standard and legacy FAPROTAX integrations and cross-OS Sankey check remain.
The final pushed SHA and its actual CI result are reported in the user handoff;
an earlier green workflow is not evidence for a later commit.

## Limits and next gate

The supplied representative FASTA/count pairs are real Illumina/DADA2 ASV
datasets, not ONT consensus datasets. Exact feature correspondence was verified,
but the native sequence trial produced no function table. After bounded native
portability fixes, APPLES changed 15 reference splits and the mapper rejected
tree matching. No compatibility check was bypassed and no non-native runtime
was introduced. The one-feature trial also approached laptop memory limits;
the complete 44-feature input was not attempted. Sequence mode has no verified
usage example. Its next gate is a successful compatible alignment/placement and
function-prediction run before exposing a feature-input production interface.

Taxonomy predictions are phenotype associations, not evidence of genes,
pathway completeness, expression, pathogenicity or activity. Their function
counts overlap. Root fallback predictions make positive-function coverage
unsuitable as a placement-accuracy estimate. Only one mapper thread, default
per-function models, no normalization and unweighted binary abundance are
verified; probability-weighted inference and general determinism are unverified.
