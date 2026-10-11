# ont-wf16S-amplicon-postprocess v0.4.9 — Astra release-validation patch

Review date: 2026-09-28. Implementation owner: Luna.
Baseline: `WhyAdr/ont-wf16S-amplicon-postprocess@843ae1ba150a807bf1c9d15ef8d3cf1f83d7f0eb`.

## Decision

**The three original taxonomy blockers are closed at the reviewed SHA.** This is a regression-coverage and release-evidence patch, not a request to reimplement the taxonomy resolver. No additional production defect in those three findings was established.

The handoff understated the release evidence: current CI already passed R testthat, validate-only no-mutation checks, release integrations and verification on Ubuntu and Windows. Rscript is missing locally, but that is a local execution limitation, not proof that the project lacks successful R validation. The annotated-tag check was skipped on this branch run and remains a tag-time requirement. No v0.4.9 tag was present in the fetched repository.

No remote changes or tags were created. The R diff below is proposed and **not locally executed**; unlike the genbank candidate, it must be validated in R before merging.

## Findings resolved against source and tests

| Concern | Verified implementation and evidence | Disposition |
| --- | --- | --- |
| Mixed unresolved outcomes mislabeled as context mismatch | `TerminalFailure("unresolved", "E_TAXONOMY_UNRESOLVED", ...)`; individual `lookup_outcomes` retain their semantics. `test_mixed_semantic_unresolved_nodes_keep_lookup_statuses_in_terminal_failure` checks not-found plus ambiguous outcomes and aggregate failure. | Closed |
| Local execution errors mislabeled as remote invalid responses | `terminal_failure_for_exception` preserves typed request/response failures and maps other exceptions to `execution_failed` / `E_TAXONOMY_EXECUTION`. Boundary tests cover OSError, SafeRequestFailure and InvalidResponse. | Closed |
| Double online pass in ordinary refresh | R seam test records one local `--validate-only` call without `--online-preflight`, followed by one ordinary execution call. Python test separately checks one entry into the online resolution phase. | Closed |
| Cleanup with spaces/metacharacters | POSIX `rmdir` argument is passed through `shQuote(..., type="sh")`; Windows uses quoted `cmd.exe /d /c rmdir /q`. Existing test covers spaces, `&`, and brackets on both OS runners. | Existing tested case closed; broaden regression below |
| Preserving evidence | Cleanup checks ownership/containment and emptiness, then uses nonrecursive removal. Existing tests retain diagnostic evidence and reject an unowned directory. | Preserve current behavior |
| Credential redaction | Python diagnostic tests and R seam assertions cover contact email, API key and credential-bearing URL redaction. | Existing contract retained |

“One online pass” means one resolver phase in ordinary execution. It does **not** mean one HTTP request: multiple taxa, ESearch/EFetch and bounded retries legitimately require multiple requests. Explicit user-requested online validation is a separate mode, not a reason to add online preflight back into ordinary execution.

Primary source links: [resolver](https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/blob/843ae1ba150a807bf1c9d15ef8d3cf1f83d7f0eb/analysis/utils/ncbi_taxonomy.py), [Python regressions](https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/blob/843ae1ba150a807bf1c9d15ef8d3cf1f83d7f0eb/tests/test_ncbi_taxonomy.py), [R invocation contract](https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/blob/843ae1ba150a807bf1c9d15ef8d3cf1f83d7f0eb/tests/testthat/test-kreport.R), [cleanup](https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/blob/843ae1ba150a807bf1c9d15ef8d3cf1f83d7f0eb/analysis/utils/atomic_io.R).

## Evidence actually obtained

- Exact repository head matched Luna's note.
- Local full Python discovery: **79 tests passed** (`python -m unittest discover -s tests -p 'test_*.py'`).
- Local Node Sankey client test: **passed**.
- Local release metadata check: **passed for v0.4.9**.
- Local `git diff --check`: passed on baseline; R was unavailable and not run locally.
- [CI run 36230903935](https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/actions/runs/36230903935): both primary Ubuntu/Windows jobs succeeded; both FAPROTAX jobs succeeded; cross-OS Sankey byte equality succeeded.
- In both primary jobs, unit/regression, all three no-mutation steps, full release integration, zero-classified cohort and phylogenetic integration steps succeeded. The workflow invokes `tests/verify_release_run.R` after each corresponding integration.
- The exact annotated-tag metadata step was **skipped**, as expected for a branch run. Do not report it as passed.

## Diff 1 — broaden cleanup regression without changing production behavior

Replace the existing single-path regression with the following loop. Its purpose is to ensure shell quoting treats apostrophes, dollar signs and backticks as literal filename characters while preserving adjacent evidence. A unique root avoids collisions between repeated test runs. Existing evidence-retention and wrong-parent tests remain unchanged.

```diff
--- a/tests/testthat/test-atomic-ownership.R
+++ b/tests/testthat/test-atomic-ownership.R
@@ -410,16 +410,26 @@
 })

 test_that("empty taxonomy diagnostics cleanup quotes paths with spaces and shell metacharacters", {
-  root <- file.path(tempdir(), "diagnostics cleanup & [shell]")
+  root <- tempfile("diagnostics cleanup ")
   dir.create(root, recursive = TRUE, showWarnings = FALSE)
   on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
-  final_root <- file.path(root, "published output")
-  diagnostics <- create_taxonomy_diagnostics_dir(final_root, new_transaction_id())
-  diagnostics_root <- dirname(diagnostics)
-
-  expect_true(cleanup_empty_taxonomy_diagnostics(diagnostics, diagnostics_root))
-  expect_false(dir.exists(diagnostics))
-  expect_true(dir.exists(diagnostics_root))
+  # These characters are legal filenames on both supported platforms.
+  # Percent expansion and delayed expansion in cmd.exe are separate concerns.
+  names <- c("spaces & [shell]", "apostrophe's path", "dollar $HOME", "backtick `literal`")
+  for (name in names) {
+    parent <- file.path(root, name)
+    dir.create(parent, recursive = TRUE)
+    sentinel <- file.path(parent, "preserve.txt")
+    writeLines("preserve sibling evidence", sentinel)
+    final_root <- file.path(parent, "published output")
+    diagnostics <- create_taxonomy_diagnostics_dir(final_root, new_transaction_id())
+    diagnostics_root <- dirname(diagnostics)
+
+    expect_true(cleanup_empty_taxonomy_diagnostics(diagnostics, diagnostics_root), info = name)
+    expect_false(dir.exists(diagnostics), info = name)
+    expect_true(dir.exists(diagnostics_root), info = name)
+    expect_identical(readLines(sentinel), "preserve sibling evidence", info = name)
+  }
 })

 test_that("taxonomy diagnostics evidence and unexpected entries are retained", {
```

This does not certify every shell/path edge case. In particular, Windows `%NAME%` expansion and delayed `!NAME!` expansion are not established safe by the current test or this proposed matrix. Do not infer universal Windows shell safety from a green spaces-and-ampersand test. If supporting such paths is required, add an isolated Windows test with a controlled environment variable and sentinel directories first, then consider an argv-safe removal primitive rather than ad hoc escaping. That is a bounded follow-up; no Windows reproducer was executed here, so it is not presented as a confirmed new release blocker.

## Diff 2 — add a release evidence checklist

Add `patch/v0.4.9-release-validation-checklist.md` using this complete new-file content (the candidate SHA and run ID are filled only after the candidate CI completes):

```diff
--- /dev/null
+++ b/patch/v0.4.9-release-validation-checklist.md
@@ -0,0 +1,31 @@
+# v0.4.9 release validation
+
+## Candidate identity
+
+- Candidate commit SHA: PENDING
+- Validation workflow run URL: PENDING
+- Reviewed baseline: 843ae1ba150a807bf1c9d15ef8d3cf1f83d7f0eb
+- Baseline run: https://github.com/WhyAdr/ont-wf16S-amplicon-postprocess/actions/runs/36230903935
+
+## Required evidence on the candidate
+
+- [ ] R 4.5.3 / Python 3.12 Ubuntu primary job succeeds.
+- [ ] R 4.5.3 / Python 3.12 Windows primary job succeeds.
+- [ ] Both FAPROTAX jobs succeed.
+- [ ] Cross-OS Sankey byte comparison succeeds.
+- [ ] Release metadata check succeeds.
+- [ ] Validate-only ordinary, Krona and Pavian checks do not create outputs.
+- [ ] Full, BAM-stats, zero-classified and phylogenetic integrations verify.
+- [ ] Cleanup path matrix and evidence-retention regressions pass.
+- [ ] Mixed unresolved outcomes retain individual statuses and aggregate unresolved.
+- [ ] Local execution, invalid response and request failures remain distinct.
+- [ ] Ordinary refresh performs one local preflight plus one execution phase.
+- [ ] Credential redaction tests pass without live credentials.
+
+## Release identity
+
+- [ ] Maintainer approves the exact validated commit for release.
+- [ ] Annotated v0.4.9 tag resolves to that commit.
+- [ ] check_release_metadata.py --require-tag v0.4.9 passes on that tag.
+- [ ] Tag-triggered validation succeeds; skipped checks are recorded explicitly.
+- [ ] Any later source change receives new candidate validation.
```

No CI rewrite is necessary: the existing workflow already runs these suites and integrations. Do not add a live NCBI dependence to mandatory CI; deterministic mocks test semantics without network outages, mutable taxonomy or credentials changing release outcomes.

## Implementation and acceptance sequence

1. Confirm the baseline SHA or reconcile newer changes. Apply the test diff and add the checklist. Leave the already-fixed resolver and ordinary refresh invocation logic intact.
2. Restore the locked R environment and run the existing commands below. If the broadened matrix fails, retain the failing path and OS and diagnose it before any cleanup implementation change.
3. Push the candidate through the existing workflow and retain the exact SHA/run association. Baseline green checks are not validation of the new test commit.
4. Complete the checklist with actual evidence. Release/tag work remains the maintainer's separate step; this handoff does not ask Luna to publish automatically.

```bash
Rscript analysis/install_packages.R --restore
Rscript tests/testthat.R
python -m unittest discover -s tests -p 'test_*.py'
node --test tests/test_taxonomy_sankey_client.mjs
python tests/check_release_metadata.py
python -m compileall -q analysis tests
git diff --check
```

Use the exact integration invocations already in `.github/workflows/ci.yml`; they include `tests/verify_release_run.R`, optional modules, synthetic fixtures and cross-OS artifact comparison. There is no need to manually repeat successful baseline integrations solely because this review environment lacks R. The new candidate needs its own normal CI run.

After the maintainer creates the annotated tag on the validated commit, the tag-time check is:

```bash
python tests/check_release_metadata.py --require-tag v0.4.9
```

Do not run that command before the tag exists and then mistake the expected failure for a taxonomy defect. Keep VERSION, CITATION.cff, README release anchor and CHANGELOG consistent; if changing the planned release date, update its paired metadata fields together.

## Handoff to Luna

Implement the regression/checklist patch, return the final SHA, CI run and any failure details, and mark the old three taxonomy issues resolved. This review gives conditional release sign-off for the reviewed behavior, subject to normal candidate/tag identity gates. It does not promise safety for arbitrary shell paths or claim local R execution.
