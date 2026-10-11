# ont-wf16S-amplicon-postprocess v0.4.9 Release-Hardening Plan

**Target repository:** WhyAdr/ont-wf16S-amplicon-postprocess
**Target release:** v0.4.9
**Audit anchor supplied with the draft:** 985b42ee352947e9bcc61846b7d1308a80202486
**Purpose:** close the remaining taxonomy-outcome, refresh-regression, and diagnostics-lifecycle gaps without redesigning the v0.4.9 refresh transaction.

## Review basis and implementation preflight

This plan has been reviewed for internal consistency against the attached audit findings and cited audit commit. Before changing code, the implementer must inspect the current main branch, confirm the audited behavior and symbols still match, and adapt the patch to any intervening commits. Treat the SHA above as the planâ€™s evidence anchor, not as permission to patch an old checkout or overwrite newer work.

The refinements below resolve the draft's open choices:

- Aggregate terminal failures are distinct from per-taxid semantic lookup outcomes.
- NCBI request and response exceptions are classified before the generic execution handler.
- Do not add an environment-controlled production test hook such as WF16S_TEST_MODE.
- The one-online-capable-pass contract is tested offline at the R orchestration and Python resolver layers.
- Diagnostics cleanup is registered inside the actual transaction function, after rollback handling is established; it is not placed at script top level.
- Cleanup applies only during ordinary R return/error unwinding. It cannot clean up after SIGKILL, power loss, or other abrupt process termination.
- The existing NCBI email/configuration contract is preserved. This patch does not move or weaken credential validation.
- Existing credential-redaction protections receive explicit regression coverage.

## 1. Summary and priorities

The supplied audit reports that v0.4.9 already has the intended refresh architecture:

- --validate-only --refresh-taxonomy performs local validation without contacting NCBI;
- --online-preflight is restricted to validation mode;
- refresh cache publication is deferred and rollback-protected;
- resolver diagnostics are written outside disposable module staging;
- Python and R redact resolver errors before surfacing them.

Do not redesign those mechanisms in this patch. Address the following items:

| Priority | Item | Required result |
|---|---|---|
| Release blocker | RH-001 | Aggregate unresolved failure reports outcome=unresolved, not context_mismatch. |
| Release blocker | RH-002 | Local resolver failures report execution_failed; actual NCBI response failures remain response_invalid. |
| Release blocker | RH-003 | Offline regression gates protect local refresh validation, online-preflight mode constraints, and one online-capable execution invocation. |
| Required regression gate | RH-004 | Empty transaction-specific diagnostic directories are cleaned safely; evidence-bearing diagnostics and rollback recovery material are retained. |
| Required regression gate | Redaction | Sentinel email/API-key/credential-bearing URL values never appear in resolver events, failure summaries, or surfaced R errors. |

## 2. RH-001: distinguish aggregate unresolved from per-lookup outcomes

The current terminal E_TAXONOMY_UNRESOLVED summary is described as using context_mismatch, even though the unresolved set may contain not_found, ambiguous, and context_mismatch results.

### Required contract

- Preserve each individual lookup's semantic result.
- Report the aggregate terminal failure as outcome=unresolved with code=E_TAXONOMY_UNRESOLVED.
- Keep SEMANTIC_OUTCOMES limited to per-lookup results. Do not add unresolved there merely to serialize a terminal summary.
- Keep existing JSON field names and required fields stable. If the current summary writer only accepts a LookupOutcome, add a narrow terminal-summary representation or argument instead of recording an aggregate as though it were a lookup result.
- Add semantic_counts only if it can be added compatibly and is useful to downstream readers. It is optional for this patch.

### Pseudocode diff

~~~diff
diff --git a/analysis/utils/ncbi_taxonomy.py b/analysis/utils/ncbi_taxonomy.py
@@
 if args.unresolved_policy == "error" and unresolved:
-    events.write_failure_summary(
-        LookupOutcome(
-            "context_mismatch",
-            code="E_TAXONOMY_UNRESOLVED",
-            message=f"{len(unresolved)} taxonomy nodes remain unresolved",
-        ),
-        len(lookup_outcomes),
-        len(lookup_outcomes),
-    )
+    events.write_failure_summary(
+        TerminalFailure(
+            outcome="unresolved",
+            code="E_TAXONOMY_UNRESOLVED",
+            message=f"{len(unresolved)} taxonomy nodes remain unresolved",
+        ),
+        processed=len(lookup_outcomes),
+        total=len(lookup_outcomes),
+    )
     return 1
~~~

Use the existing summary API if it can represent a terminal outcome without adding unresolved to per-lookup classification. TerminalFailure above is illustrative, not a required class name.

### Regression

Construct a case with a mixed unresolved set (at least two different semantic results). Assert:

- each event retains its original per-lookup status;
- the terminal JSON summary has code=E_TAXONOMY_UNRESOLVED and outcome=unresolved;
- no terminal summary claims the whole set was context_mismatch.

## 3. RH-002: classify response failures before local execution failures

The draft's generic exception change is incomplete unless exception precedence is explicit. A malformed NCBI response must not fall through to the generic handler merely because its parser exception is a ValueError or OSError.

### Required classification

| Failure source | Code | Terminal outcome |
|---|---|---|
| NCBI transport/request failure | E_NCBI_REQUEST | request_failed |
| Malformed, inconsistent, or explicit-error NCBI response | E_NCBI_RESPONSE | response_invalid |
| Aggregate policy failure after resolution | E_TAXONOMY_UNRESOLVED | unresolved |
| Local config, input, filesystem, or resolver execution failure | E_TAXONOMY_EXECUTION | execution_failed |

Catch the project's typed request and response exceptions before local/generic exceptions. If InvalidResponse is not already the exclusive representation of response-validation failures, normalize those raise sites to that type in this patch. Do not reclassify HTTP/transport errors as execution failures.

At the command boundary, expected local failures and otherwise-unhandled ordinary exceptions may be reported as execution_failed so the CLI exits nonzero with a structured summary. Do not catch BaseException; preserve normal handling of interrupts and process exits. Keep surfaced messages redacted and do not add raw traceback or response bodies to user-facing diagnostics.

### Pseudocode diff

~~~diff
diff --git a/analysis/utils/ncbi_taxonomy.py b/analysis/utils/ncbi_taxonomy.py
@@
+except NCBIRequestError as exc:
+    outcome = failure("request_failed", "E_NCBI_REQUEST", exc)
+    write_safe_failure(outcome)
+    return 1
+except InvalidResponse as exc:
+    outcome = failure("response_invalid", "E_NCBI_RESPONSE", exc)
+    write_safe_failure(outcome)
+    return 1
 except (OSError, ValueError, json.JSONDecodeError) as exc:
-    outcome = failure("response_invalid", "E_TAXONOMY_EXECUTION", exc)
+    outcome = failure("execution_failed", "E_TAXONOMY_EXECUTION", exc)
     write_safe_failure(outcome)
     return 1
+except Exception as exc:
+    # Only at the CLI boundary; retain a safe, structured fatal result.
+    outcome = failure("execution_failed", "E_TAXONOMY_EXECUTION", exc)
+    write_safe_failure(outcome)
+    return 1
~~~

Adapt exception names and catches to the actual implementation. Avoid redundant json.JSONDecodeError in the tuple if it is already covered by ValueError; the important behavior is explicit NCBI exception precedence.

### Regression set

1. A local input/hash/configuration failure yields E_TAXONOMY_EXECUTION / execution_failed.
2. A malformed NCBI response remains E_NCBI_RESPONSE / response_invalid.
3. A request/transport failure remains E_NCBI_REQUEST / request_failed.
4. All fatal cases return nonzero and produce summaries with the existing schema.
5. R summary rendering preserves the code and safe message for execution_failed, unresolved, response_invalid, and request_failed.

## 4. RH-003: protect refresh-vs-validation semantics offline

Release CI must not depend on live NCBI. Use a fixture with a real source taxonomy cache that lacks at least one assignment from the test abundance table. Use a unique sentinel email in the test environment and keep the fixture/cache in a temporary directory.

### 4.1 Local refresh validation

Run the real R controller process with:

~~~text
--config <fixture-config>
--modules kreport
--refresh-taxonomy
--validate-only
--output-dir <temporary-output>
~~~

Test both unresolved_policy: warn and unresolved_policy: error. In both cases local validation must exit zero and report taxonomy_resolution=pending_online; the unresolved policy is enforced only after an online-capable resolution attempt.

Assert before and after the invocation:

- zero NCBI requests;
- source cache bytes are identical;
- no taxonomy journal or source-cache backup was created;
- requested output directory does not exist;
- .wf16s-diagnostics was not created for this validation;
- the output indicates validation passed and pending_online.

The test must exercise the R controller, not only invoke ncbi_taxonomy.py directly. Keep the current NCBI email/configuration requirement unchanged; supply a valid sentinel email if required.

### Pseudocode

~~~diff
diff --git a/tests/testthat/test-transaction-process.R b/tests/testthat/test-transaction-process.R
@@
+test_that("refresh validate-only stays local and leaves no outputs or cache mutations", {
+  fixture <- create_unresolved_refresh_fixture(unresolved_policy = "error")
+  before <- fingerprint_cache_and_transaction_sidecars(fixture$cache)
+
+  result <- run_transaction_process(c(
+    "--config", fixture$config,
+    "--modules", "kreport",
+    "--refresh-taxonomy",
+    "--validate-only",
+    "--output-dir", fixture$output
+  ), env = c(NCBI_EMAIL = "redaction-test@example.invalid"))
+
+  expect_identical(result$status, 0L)
+  expect_match(paste(result$stdout, result$stderr), "pending_online")
+  expect_match(paste(result$stdout, result$stderr), "Validation check PASSED")
+  expect_identical(fingerprint_cache_and_transaction_sidecars(fixture$cache), before)
+  expect_false(dir.exists(fixture$output))
+  expect_false(dir.exists(fixture$diagnostics_root))
+})
~~~

Use existing test helpers and path conventions. Include the output parent/root in the fixture so the diagnostic-directory assertion checks the actual location the controller would use.

### 4.2 Online-preflight mode matrix

Protect the full mode contract:

| Flags | Expected |
|---|---|
| --refresh-taxonomy --validate-only | Local-only validation; no request. |
| --refresh-taxonomy --validate-only --online-preflight | Accepted online validation; controlled resolver stub; no cache/output publication. |
| --refresh-taxonomy --online-preflight | Rejected before output/cache mutation. |
| --validate-only --online-preflight without refresh | Rejected before output/cache mutation. |

For the accepted online-preflight case, use a deterministic fake resolver or mocked Python request path. Assert that a fatal request/response result propagates nonzero and does not publish a candidate cache or output tree. Do not call live NCBI from CI. Preserve any diagnostic evidence that the online-preflight attempt intentionally creates; do not assert that its diagnostic path is absent unless that is the current documented contract.

The validation-mode error identifier should be asserted exactly where stable (for example, E_ONLINE_PREFLIGHT_MODE); otherwise assert the user-facing mode error and nonzero exit.

### 4.3 One online-capable execution invocation, without a production test hook

Do not add an environment-variable-controlled test probe to production code. Test the contract in two deterministic layers:

1. **R orchestration test:** stub/inject the subprocess runner at the narrow resolver-invocation boundary. For ordinary refresh execution, assert exactly two resolver invocations:
   - local preflight: refresh + validate-only, without online-preflight;
   - execution: refresh execution, without online-preflight.

   This means one online-capable execution invocation, not one HTTP request total. The number of NCBI HTTP calls can vary with unresolved taxa and batching.

2. **Python resolver test:** with network calls mocked, prove validate-only refresh makes zero requests and the execution path enters the online resolver phase once. Assert the expected per-taxon/batch calls separately from the number of resolver phases.

If there is no runner seam, introduce a small resolver invocation helper with an injected runner argument defaulting to the existing subprocess implementation. Do not add WF16S_TEST_MODE, a test log environment variable, or a user-facing config option solely for tests. Keep the R process regression in Â§4.1 and the CLI mode matrix in Â§4.2 even if Â§4.3 is implemented as an orchestration-level test rather than a full pipeline run.

### Pseudocode for the orchestration assertion

~~~text
calls = run_refresh_controller(fixture, runner = recording_stub)

assert calls == [
  {refresh: true, validate_only: true, online_preflight: false},
  {refresh: true, validate_only: false, online_preflight: false}
]
~~~

## 5. RH-004: safe diagnostics cleanup and rollback lifecycle

Diagnostics are written outside disposable module staging so that resolver evidence survives stage cleanup and output rollback. Keep that property.

### Cleanup contract

- Create the normal transaction-specific diagnostics directory for execution-mode refresh.
- Local-only --validate-only --refresh-taxonomy must not create that directory or its shared root.
- Explicit online-preflight may need isolated diagnostic output for a resolver failure. Preserve the current tested behavior, but do not create the normal execution transaction/staging tree or publish cache/output. Document and test that path separately rather than treating it as local-only validation.
- Delete only an empty transaction-specific directory owned by the current transaction.
- Never recursively delete a diagnostics directory.
- Preserve any non-empty diagnostics directory, including taxonomy_events.jsonl, taxonomy_failure.json, and unexpected files.
- Never remove the shared .wf16s-diagnostics root in this patch.
- Run cleanup only after the transaction's rollback/finalization logic has completed.
- Register cleanup inside the actual transaction/controller function. Do not place on.exit() at script top level.
- If errors use normal R stack unwinding, on.exit() or tryCatch(..., finally=...) may run cleanup. If the controller calls quit() directly, place cleanup in the structured transaction wrapper or return through that wrapper first.
- Cleanup cannot run after SIGKILL, power loss, or other abrupt termination. A stale empty transaction directory after such a crash is acceptable; do not add age-based pruning here.

### Pseudocode helper

~~~r
cleanup_empty_taxonomy_diagnostics <- function(path, expected_parent) {
  if (is.null(path) || !dir.exists(path)) return(invisible(FALSE))

  # Require a transaction-owned child of the expected diagnostics root.
  if (!identical(normalizePath(dirname(path), mustWork = TRUE),
                 normalizePath(expected_parent, mustWork = TRUE))) {
    return(invisible(FALSE))
  }

  entries <- list.files(
    path,
    all.files = TRUE,
    no.. = TRUE,
    recursive = FALSE,
    include.dirs = TRUE
  )
  if (length(entries) != 0L) return(invisible(FALSE))

  status <- unlink(path, recursive = FALSE, force = TRUE)
  invisible(status == 0L && !dir.exists(path))
}
~~~

Adapt parent/path validation to the project's existing transaction ownership model. The helper should fail closed: on uncertainty or unlink failure, retain the directory and emit a safe warning through the existing diagnostics/error channel.

Register it inside the function that owns the transaction and only after rollback handling has been arranged:

~~~r
run_pipeline_transaction <- function(...) {
  diagnostics_dir <- create_taxonomy_diagnostics_dir(...)
  on.exit(
    tryCatch(
      cleanup_empty_taxonomy_diagnostics(diagnostics_dir, diagnostics_root),
      error = function(e) warning("Could not clean empty taxonomy diagnostics directory",
                                  call. = FALSE)
    ),
    add = TRUE
  )

  # Existing staged execution, rollback, and transaction closure.
}
~~~

The pseudocode is illustrative: use the actual transaction function and ensure its cleanup runs after rollback. Do not register this at top level in analysis/00_run_pipeline.R.

### Diagnostics and rollback regression tests

Add or extend tests to prove:

1. An empty transaction directory is removed on ordinary early return/error before resolver evidence is written.
2. A directory containing a resolver event or failure summary is preserved.
3. A directory containing any other entry is also preserved.
4. Resolver evidence survives module-stage cleanup and downstream failure.
5. If a refresh candidate has been prepared and a later module fails, the source taxonomy cache is restored byte-for-byte and diagnostics remain available.
6. Existing rollback ownership protections still refuse to overwrite a cache changed by another actor; recovery material/journals remain when rollback cannot safely complete.
7. Successful closure removes transaction journals as before; successful diagnostic retention follows the existing documented policy.

Do not change cache publication, rollback ordering, or journal semantics unless one of these tests exposes a defect. Any such defect is a separate implementation finding and must be reported before broadening scope.

## 6. Credential redaction regression gate

The attached audit says redaction is already implemented in Python and at the R boundary. This patch should preserve it and make it testable.

Use distinct sentinels for:

- the configured contact email;
- an NCBI API-key query parameter;
- URL user-info or another credential-bearing URL form currently handled by the redactor.

Drive a controlled failure whose exception text contains those sentinels. Assert none occur in:

- taxonomy_events.jsonl;
- taxonomy_failure.json;
- Python stderr/CLI error output;
- the error surfaced by the R controller.

Also assert that the error remains useful: it reports the correct failure code/outcome and a safe explanation. Extend existing redaction tests rather than duplicating fixtures. Do not log raw exception text first and redact only the final summary.

## 7. CI and documentation

### CI

Prefer the testthat process tests for the refresh-specific validation gate because they already run under the Ubuntu/Windows matrix and can assert filesystem state portably. Avoid a second shell-only fixture step with placeholder paths and grep if the same assertions are covered in testthat.

If a separate named CI step is useful for visibility, have it invoke the real checked-in test helper and use repository-supported shell/path conventions. It must use an offline fixture and must not run --online-preflight against live NCBI.

Keep the existing Linux and Windows CI, FAPROTAX, and cross-OS Sankey jobs. Confirm the new tests run on both platforms.

### README

Document only behavior verified by tests:

- local --validate-only --refresh-taxonomy makes no NCBI request and reports pending_online when unresolved taxa remain;
- unresolved_policy: error is deferred until online resolution;
- non-empty transaction diagnostics are retained outside disposable staging;
- empty transaction-specific directories may be removed during normal cleanup;
- abrupt process termination may leave an empty transaction directory.

Do not change the current NCBI email requirement as part of this patch. Document whether local validation requires the configured email only if that behavior is verified from the implementation and test.

### CHANGELOG

Under v0.4.9, summarize:

- aggregate unresolved summaries now use unresolved;
- local failures now use execution_failed, while NCBI request/response classifications remain distinct;
- offline refresh-validation/orchestration/diagnostics regression coverage;
- safe empty diagnostics cleanup and retention of evidence.

Do not bump beyond 0.4.9; do not tag until every release gate passes.

## 8. Required test and release matrix

Use the repository's canonical test commands/CI scripts. At minimum, run the equivalent of:

~~~bash
python -m compileall -q analysis tests
python -m unittest -v tests/test_ncbi_taxonomy.py tests/test_check_committed_whitespace.py
Rscript tests/testthat.R
python tests/check_release_metadata.py
git diff --check
~~~

Then run targeted offline checks for:

- local refresh validation with unresolved_policy: warn;
- local refresh validation with unresolved_policy: error;
- invalid online-preflight flag combinations;
- accepted online-preflight with a controlled resolver failure;
- R orchestration call sequence for normal refresh;
- Python local-vs-online resolver behavior with requests mocked;
- redaction sentinels across Python and R surfaces;
- diagnostics cleanup, evidence retention, cache rollback, and rollback conflict protection.

Run existing release integration checks that are feasible and already part of the projectâ€™s release gate, including the existing output verifier. Do not make the release gate depend on live NCBI.

Record exact test counts, skips, warnings, and final commit SHA. Confirm Ubuntu and Windows CI results on that exact SHA, plus FAPROTAX and cross-OS reproducibility jobs.

## 9. Acceptance criteria

### Outcome classification

- E_NCBI_REQUEST summaries use request_failed.
- E_NCBI_RESPONSE summaries use response_invalid.
- E_TAXONOMY_UNRESOLVED summaries use unresolved.
- E_TAXONOMY_EXECUTION summaries use execution_failed.
- Per-lookup not_found, ambiguous, and context_mismatch semantics are unchanged.
- unresolved_policy=warn continues to tolerate semantic unresolved results.
- unresolved_policy=error is applied only after an online-capable resolution attempt in execution/online-preflight modes.

### Refresh semantics

- Ordinary refresh preflight is local and performs no NCBI request.
- Ordinary refresh execution does not pass --online-preflight.
- Only the combination --validate-only --refresh-taxonomy --online-preflight enables explicit online validation.
- Invalid --online-preflight combinations fail before output/cache publication.
- Local refresh validation exits successfully with pending_online for unresolved taxa, even with unresolved_policy=error.
- Local validation does not mutate the source cache or create output/diagnostic transaction directories.

### Transaction and diagnostics safety

- Resolver failure does not publish or mutate the source taxonomy cache.
- Deferred candidate commit and rollback behavior remains intact.
- Later-module failure restores the original cache when ownership checks permit.
- Rollback refuses to overwrite cache changes from another actor.
- Diagnostic evidence survives module-stage cleanup and rollback.
- Only empty transaction-specific diagnostics directories are eligible for cleanup.
- Cleanup never deletes the shared diagnostics root or non-empty evidence.
- Journals/backups/recovery material follow existing transaction contracts.

### Redaction and release gates

- Sentinel email, API key, and credential-bearing URL values are absent from diagnostic files and surfaced errors.
- Linux and Windows CI pass on the final audited commit.
- FAPROTAX and cross-OS reproducibility jobs remain green.
- Python taxonomy tests and full testthat suite pass with no new untriaged warnings.
- Release metadata checks report 0.4.9.
- git diff --check passes.
- No tag is created before final release audit.

## 10. Explicit non-goals

Do not expand this patch into:

- asynchronous taxonomy refresh or parallel NCBI requests;
- persistent/distributed NCBI caching or Entrez redesign;
- taxonomy-cache schema migration or rank-system redesign;
- manifest schema revision;
- age-based diagnostics pruning or crash-recovery scavenging;
- Krona, Pavian, Sankey, FAPROTAX, or abundance-accounting changes;
- a production environment/configuration switch solely to facilitate tests;
- any new network dependency in CI.

## 11. Handoff instructions

Implement against the current main branch after confirming the audit anchor and current call sites.

1. Fix terminal taxonomy outcome representation and exception precedence.
2. Add focused Python outcome, request/response, and redaction regressions.
3. Add R process tests for local refresh validation and the full online-preflight mode matrix.
4. Add an offline R orchestration test for the one online-capable execution invocation; use runner injection, not a production test-mode environment variable.
5. Add safe empty-directory cleanup inside the transaction function and test evidence retention plus rollback behavior.
6. Update README and CHANGELOG only for verified behavior.
7. Run the full cross-platform test/release matrix.
8. Inspect the final diff for unrelated changes and verify release metadata remains 0.4.9.
9. Report final commit SHA, exact test counts, skips/warnings, Ubuntu/Windows and other CI results, and any deferred item.

Do not tag v0.4.9 until the final audited commit satisfies every acceptance criterion. Hand that exact HEAD back for one last release audit before tagging.
