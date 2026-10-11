# wf-16s Kraken2/Bracken support: reviewed implementation handoff

Reviewed: 2026-10-10. Proposed target: v0.5.0.
Verified local baseline: `3522ebce073ed44573183e2a88da833da9c6636b`, VERSION `0.4.9`.
Workspace: `D:/W/AAy_Amplicon`.
Original draft: `kraken2-support-implement.md` (preserved).
Status: ready for implementation against the contracts below; implementation and release validation remain pending.

## 1. Decision and corrections to the original draft

Add a classifier-specific input adapter and accounting model. Expanding accepted field counts alone will not make these inputs valid: Bracken estimates are not the same population as Kraken read assignments.

| Draft assumption | Audit finding and required correction |
| --- | --- |
| Biomecon abundance totals equal 40,204 assignments and 15,559 classified reads | False. Its abundance table sums to 38,629, including 13,984 species estimates and 24,645 raw unclassified reads. |
| Species-matrix semantics are identical across classifiers | Matrix arithmetic can be reused, but Kraken species rows contain Bracken estimated counts. Read-accounting denominators, labels, provenance, and verification must change. |
| Detect classifier by accepting either 5 or 6 fields per row | Select the schema from validated `params.classifier`. Reject mixed formats and classifier/schema disagreement throughout a file, including across chunk boundaries. |
| Updating one R validation chunk and Python's lineage index is enough | R scans and then assembles the file in a second pass. Both passes, all callers, Python validation, QC assertions, manifests, and independent release verification must agree. |
| Both bundle-level params files are suitable contracts | Use the inner `output/params.json`. The launch-level `params.json` lacks `database_sets`, `bracken_length`, and other defaults. Do not synthesize missing provenance. |
| Bracken `.report` is a six-column Kraken report | The supplied Bracken files are headered seven-column abundance tables. Dispatch by the documented producer schema, never by extension alone. |
| Species names can join the sources | Some Bracken names are older names for the same TaxIDs. Use verified TaxID-to-lineage mappings; preserve alias evidence. |
| Two samples should yield 2D PCoA and NMDS figures | Preserve current rank/sample gates. One-dimensional scores and explicit skips are correct for this pair. |
| ITS data should run FAPROTAX unchanged as an ordinary functional analysis | All positive abundance rows here are Eukaryota. The existing no-eligible-prokaryote skip is the correct outcome when FAPROTAX is requested. |

This handoff authorizes no release publication by itself. Implement and validate locally; commit, push, tag, or publish only under the execution chat's authorization. Keep the supplied bundles and unrelated plans out of staging. Never upload sequences or invoke remote BLAST.

## 2. Measured input evidence

These are full-file measurements, not estimates from file previews. Reproduce them with:

```powershell
python patch/kraken2-input-audit.py
```

The script reads metadata, assignments, abundance, and reports locally, sends JSON to stdout, and makes no network requests. `patch/kraken2-input-audit.json` records source paths, SHA-256 hashes, aliases, and placeholder differences. It is audit evidence, not a production parser; its local placeholder comparison must not be copied into a general resolver without the stricter rules in section 6.

| Quantity | Biomecon_ITS-kraken2 | Knegative_ITS-kraken2 |
| --- | ---: | ---: |
| Assignment rows, N | 40,204 | 26,643 |
| Kraken classified / positive TaxID, C | 15,559 | 5,469 |
| Kraken unclassified / TaxID 0, U | 24,645 | 21,174 |
| C + TaxID 0 | 0 | 0 |
| Bracken species estimate sum, E | 13,984 | 4,796 |
| Abundance unclassified | 24,645 | 21,174 |
| Abundance table sum, T = E + U | 38,629 | 25,970 |
| Accounting difference, D = C - E = N - T | 1,575 | 673 |
| Positive species rows | 68 | 35 |
| Assignment read lengths | 250–1,000 bp | 250–996 bp |

All assignment rows have exactly six fields and unique read IDs. The raw Kraken report root, unclassified, and per-TaxID direct counts agree with assignments in both bundles. Every Bracken row's `new_est_reads` agrees with its abundance row after a unique TaxID-backed lineage match. All abundance lineages have eight ranks; positive assignment lineages have seven, omitting kingdom.

Both complete output contracts specify: `classifier=kraken2`, `database_set=ncbi_16s_18s_28s_ITS`, `taxonomic_rank=S`, `min_len=250`, `max_len=1000`, `min_read_qual=15`, `kraken2_confidence=0`, `bracken_threshold=10`, `bracken_length=null`, `abundance_threshold=1`, and assignments enabled. Their recorded identity/coverage values are 95/90, but those values are inactive for this classifier. Tool metadata records Kraken 2.1.3 and Bracken v3.0.1. Local stdout records wf-16s v1.6.1 and short revision `e2b77380f8`.

The accounting difference is observed, but this audit does not identify how much comes from threshold exclusion, redistribution limits, or integer estimation. Do not call D failed alignment reads, newly unclassified reads, contaminants, or individually identified omitted reads. Do not infer the effective Bracken model length from a null request.

### Concrete local paths

For each sample, the producer root is nested twice:

```text
wf-16s_Biomecon-ITS-kraken2/wf-16s_Biomecon-ITS-kraken2/output/
wf-16s_Knegative-ITS-kraken2/wf-16s_Knegative-ITS-kraken2/output/
```

Use the respective files beneath that root:

```text
params.json
abundance_table_species.tsv
reads_assignments/<SampleID>_lineages.kraken2.assignments.tsv
kraken2/<SampleID>.kraken2.report.txt
bracken/<SampleID>.kraken2_bracken.report
versions_all.txt
```

### Primary-source boundaries

Native Kraken output has five fields; wf-16s adds lineage as field six. Numeric hit runs and `A:<count>` must be recognized; paired and protein formats use additional separators outside this initial ONT single-read scope. See the [Kraken manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format).

Bracken estimates abundance at a chosen rank and applies threshold-dependent redistribution. Its optional Kraken-style report is a different artifact from its headered abundance output. See the [Bracken documentation](https://github.com/jenniferlu717/Bracken#step-3-run-bracken-for-abundance-estimation).

The [wf-16s v1.6.1 entry point](https://github.com/epi2me-labs/wf-16s/blob/v1.6.1/main.nf) delegates Kraken work to wf-metagenomics. That tag links [submodule commit b0607906fa105301599f8c1f0520dd5db5ed3dc4](https://github.com/epi2me-labs/wf-metagenomics/tree/b0607906fa105301599f8c1f0520dd5db5ed3dc4). Its [database preparation source](https://github.com/epi2me-labs/wf-metagenomics/blob/b0607906fa105301599f8c1f0520dd5db5ed3dc4/modules/local/databases.nf) can choose a distribution length from database contents when the requested model is unavailable.

The local short revision could not be resolved through the public fetches in this audit; the detailed Kraken subworkflow and abundance-building source were also unavailable through those fetches. The public tag is supporting evidence, not a verified exact checkout of these runs. Before generalizing beyond these fixtures, Gemini must inspect the matching producer/submodule source, record exact revisions or report the limitation, and confirm estimate rounding, rank filling, and table generation. Never substitute current upstream code silently. Observed file contracts and measured counts above remain independently verified.

## 3. Supported scope and parameter contract

Initial support covers wf-16s single-read, DNA, lineage-enriched Kraken2 assignments plus species-level Bracken abundance using the two maintained NCBI database sets. Keep existing rejections for SILVA/Greengenes namespaces, custom reference/taxonomy/database overrides, unsupported ranks, and unknown classifiers. Native Kraken files without added lineage, paired reads, translated protein hit lists, and arbitrary Bracken-only producers remain outside scope.

Retain common producer-field validation. For Kraken additionally require finite `kraken2_confidence` in [0,1], nonnegative safe whole-number `bracken_threshold`, and a present `bracken_length` key whose value is null or a positive safe integer. Reject missing, malformed, Boolean, array, and nonfinite values. For minimap2, retain its existing required identity/coverage contract without requiring Kraken fields. Validate recorded inactive values if present, and expose their inactive status explicitly.

Extend `extract_upstream_contract()` to record requested Kraken/Bracken settings, `abundance_source=bracken_species_estimates`, the assignment format, and threshold applicability. Distinguish requested versus effective Bracken length; record effective length as unknown unless a fingerprinted producer artifact actually supplies it. EPI2ME agent version is not workflow version. Optional version/log evidence must be fingerprinted if consumed; leave workflow version/revision null when evidence is not configured.

## 4. Input API and accounting model

Add documented optional config mappings, resolved relative to the config file like existing inputs:

```yaml
input:
  # Existing abundance_table, params_json, assignments remain in place.
  kraken2_reports: null  # SampleID -> original six-column Kraken report
  bracken_reports: null  # SampleID -> headered seven-column Bracken abundance
  params_by_sample: null  # SampleID -> complete original output/params.json
```

For Kraken runs, both report mappings are REQUIRED for every selected sample, even when QC is not requested. Fail before output staging on missing, extra, duplicate, unsafe, or classifier-inconsistent mappings. Do not infer reports through an unconstrained recursive search. Retain optional assignments: without assignments, reports still support population accounting, while read-level diagnostics are unavailable. If row identity cannot be verified without assignments, require a verified taxonomy cache/mapping rather than falling back to an unchecked name join. Reject ambiguous or unresolved count-to-row mappings independently of kreport's warning policy.

For minimap2, the new report mappings default to null and do not change behavior. Reject configured Kraken/Bracken reports on a minimap2 contract. `params_by_sample` is an additive cohort provenance option; the existing `params_json` remains the representative validated contract. If supplied, require exact sample coverage and compare normalized analysis-relevant contracts against the representative, allowing run-specific paths/IDs to differ. This permits a safely merged bundle without forging a single upstream run.

Create a small classifier adapter, preferably `analysis/utils/classifier_contract.R`, with pure functions for report parsing, contract comparison, and per-sample accounting. Keep classifier branching concentrated there and in input ingestion, QC, and provenance. Do not scatter `is_kraken` checks into every numerical algorithm.

### Minimap2: preserve current equalities

For configured assignments, retain exact total/classified/unclassified reconciliation, `effective_classified = taxid != "0"`, the allowed C+TaxID0 state, optional bamstats partitioning, and existing threshold checks. Do not weaken this branch to accommodate Kraken.

### Kraken2/Bracken: validate separate equalities

Use distinct quantities:

```text
N = number of Kraken input reads after upstream preprocessing
C = original Kraken classified reads (all positive TaxIDs, any rank)
U = original Kraken unclassified reads
E = sum of retained species new_est_reads from the Bracken abundance output
T = sum of the supplied abundance column, including its synthetic U row
D = C - E, a count difference between populations/estimates

N = C + U
assignment rows = N                         [when assignments are supplied]
assignment status-C count = C               [when assignments are supplied]
assignment status-U count = U               [when assignments are supplied]
abundance unclassified = U
abundance classified species sum = E
T = E + U
N - T = D
```

The raw report's root clade count supplies C and its unclassified row supplies U. Validate all six columns, canonical TaxIDs, safe integral counts, unique taxon rows, rank codes, indentation/tree relationships, direct-count/clade-count consistency, and rounded percentages against N. Do not sum every clade count: ancestor/descendant counts overlap. Compare per-TaxID direct counts to assignments when supplied, in addition to root totals. Support legitimate absent zero-count U/root rows only through an explicit zero-case contract and fixtures; never interpret a missing positive population as zero.

The headered Bracken schema is exactly `name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads`. Validate unique positive TaxIDs, species rank, safe counts, and row identities. Fractions are rounded estimates from Bracken and must not be required to sum to exactly 1 or equal fractions of the exported integer matrix without producer-source evidence. Verify `new_est_reads` against every corresponding abundance row, including zero/missing-row rules established from producer source. Reject extra/missing positive species and altered counts. The current files also satisfy `kraken_assigned_reads + added_reads = new_est_reads`; confirm the versioned contract before enforcing this generally.

For this supported adapter, require 0 <= E <= C; unexplained excess is an unsupported producer/accounting result, not something to clip or normalize. Test E=C and E<C. If source inspection establishes a legitimate excess from a different producer/version, add a separately reviewed contract and fixtures rather than relaxing the check globally.

Do NOT disable reconciliation, compare only U, invent a tolerance, add D to the unclassified row, rescale estimates to C, create fake read assignments, or silently replace the user's species matrix with raw Kraken species counts.

### Context and exported count semantics

Retain the supplied matrix unchanged. Existing `sample_stats$TotalReads`/`ClassifiedReads` may keep their matrix-derived numeric values for compatibility, but for Kraken they denote T/E, not N/C. Add an explicit count basis and a separate `sample_accounting` structure carrying N/C/U/E/T/D and evidence availability. Ensure validate-only summaries expose both raw reads and matrix totals. Document legacy names wherever they are exported.

Diversity, rarefaction, composition, distances, and shared taxa can retain their mathematics after regression checks. Their Kraken results describe the estimated species-count matrix. Rarefaction of estimates is a resampling analysis of that matrix, not uncertainty propagation through Bracken or sampling of observed species-resolved reads. Do not present classifier agreement or ecological significance as established by successful execution.

## 5. Assignment readers: one explicit classifier contract

Pass `classifier` through `validate_assignment_chunk()`, `scan_assignments_file()`, `read_assignments_file()`, `build_context()`, and the QC fallback call. Keep the default minimap2 where compatibility with existing direct helper calls requires it; pipeline calls must pass the validated value explicitly.

| Field | Minimap2 enriched TSV | Kraken2 enriched TSV |
| --- | --- | --- |
| Required physical fields | exactly 5 | exactly 6 |
| 1–3 | status, read ID, TaxID | status, read ID, TaxID |
| 4 | existing supported minimap2 length encodings | positive plain integer length within existing parser bounds |
| 5 | lineage | space-delimited hit runs |
| 6 | absent | lineage |
| Status/TaxID | U implies 0; C may have 0 after QC | U iff TaxID 0; C requires positive TaxID |

Keep physical line numbers, blank-line handling, trailing-empty-field preservation, compressed input support, duplicate-ID checks across chunks, canonical TaxID bounds, control-character rejection, positive-lineage validation, and before/after fingerprints across both passes. Validate homogeneous field counts before `rbind`; the draft's variable-width `rbind` can recycle shorter rows. Its `ncol(fields)` test also operates on a list before conversion. Replace those snippets with a proper parser rather than copying them.

For Kraken hit lists, validate token grammar, not just TSV-safe text. Accept numeric TaxID/count runs and `A:<count>` with bounded nonnegative integers and the supported space separator. Reject malformed runs, tabs/control characters, boundary whitespace, `|:|`, `-:-`, `--use-names` TaxID text, and paired-length encodings with a clear unsupported-format error. Do not treat hit-list TaxIDs as assigned read TaxIDs, infer lineage from hits, or require every k-mer to identify a species. Retain field six physically even when an unclassified lineage is empty. Do not retain large hit lists in the QC frame unless actually needed.

Extend Python `read_assignment_taxids()` with the same explicit classifier contract and `--classifier` CLI argument. Pass it from BOTH `run_module_preflight()` and `run_kreport()`. Preserve expected-input hash checks, conflict diagnostics, canonical-ID/length checks, and resolver failure semantics. Existing minimap2 direct calls and CLI invocations remain compatible. Add cross-language valid/invalid fixture parity; changing only Python's lineage index is insufficient.

## 6. Taxonomy identity and placeholders

Do not match the Bracken scientific-name column to abundance species text by equality. Verified examples:

| TaxID | Bracken name | Assignment/abundance species name |
| --- | --- | --- |
| 2499453 | Stomiopeltis sinensis | Exopassalora sinensis |
| 52256 | [Candida] ethanolica | Pichia ethanolica |
| 1206149 | Galactomyces pseudocandidus | Geotrichum pseudocandidum |

There are also source-specific missing-rank renderings: assignment strings such as `unclassified Dipodascales family` correspond to abundance `Dipodascales_Incertae_sedis`. The audit records Biomecon's [Candida] tunisiensis (TaxID 1183333) and Polychytrium aggregatum (110093), and Knegative's Paramicrosporidium saccamoebae (1246581). See the complete paths in the JSON evidence.

Build a unique TaxID-to-lineage map from validated assignments and/or fingerprinted taxonomy evidence. Match the eight-rank table against the seven-rank producer lineage with the existing kingdom omission plus narrowly defined missing-rank normalization confirmed against producer source. Preserve rank positions and all named ancestors. Treat inferred/filled placeholder ranks as missing nodes, never as independently resolved taxa; exclude them from inappropriate NCBI lookups. Report the original and normalized paths and their evidence. Reject named-ancestor conflicts, ambiguous TaxID mappings, and same-name/different-lineage collisions. Do not globally strip `unclassified` text or `_Incertae_sedis` suffixes to make rows pass.

The existing resolver seeds species TaxIDs by exact normalized assignment lineage; that path needs regression coverage for these placeholder differences as well as aliases. Bracken species absent from direct assignments are legitimate in principle: use verified cache/taxonomy evidence or the existing explicit taxonomy-resolution route, never invent a TaxID or suppress a count-integrity failure under `unresolved_policy=warn`.

Prepare a dedicated ITS taxonomy cache for real-data kreport validation. Do not reuse the AAy bacterial fixture cache as if it covers these fungi. Defaults stay `cache_only`; cache resolution gaps must be reported. An execution explicitly requesting `--refresh-taxonomy` may use the established NCBI name/TaxID route and transactional cache protections; no sequence upload is needed. This review did not perform a taxonomy refresh. Required strict accounting/row-identity resolution must complete before output staging, including during validate-only.

## 7. QC, exports, and downstream behavior

QC must distinguish raw Kraken read outcomes from Bracken estimates in all tables and captions, not only Figure 1c:

- Raw accounting plots use N with C/U; the center label cannot use T. Read-length strata use raw assignment classification, never Bracken species estimates.
- Add classifier, count basis, raw N/C/U, E/T/D, availability, and separate read/report and abundance/Bracken verification fields to accounting/reconciliation outputs. Define `ReconciliationPass` as all applicable contract checks passing, not N=T. An absent assignment comparison is unavailable, not a fabricated pass; report checks can still pass independently.
- C+TaxID0 is zero only for a successfully parsed Kraken assignment file; with assignments absent it is unavailable. Alignment failure partitions and identity/coverage diagnostic values are NA/not applicable for Kraken.
- Guard bamstats discovery in `build_context()`, partitioning in preflight, AND partitioning in `run_qc()`. A Kraken run must not consume stray bamstats found under a broad root. `BamstatsAvailable=FALSE` for the adapter; report alignment diagnostics as not applicable.
- Replace final QC `stopifnot` equalities that currently require RawC+RawU=AbundanceTotal and TaxID_GT0=AbundanceClassified with the explicit classifier model. Keep the existing minimap2 assertions intact.
- Use “Kraken classified” and “Kraken unclassified”; remove “Never aligned” and “aligned + passed QC” labels from Kraken plots. Display D separately as a count difference, without assigning it to individual reads.
- Set both real ITS configs' QC display range to 250–1,000 bp. The current 1,200–1,800 default would hide these reads. Replace the hardcoded “Full-length 16S” violin subtitle with wording suitable for the configured input.

Keep module 07's generated `.kreport` derived from the supplied abundance matrix. For Kraken it is a **Bracken-derived species-estimate tree** with root E, unclassified U, and percentage denominator T. It will not reproduce the original Kraken root C or total N. Preserve the original report as distinct fingerprinted input; never overwrite it. Carry this count basis through kreport, Krona, built-in explorer, Sankey JSON/HTML captions and provenance. Retain official Pavian six-column interoperability and the existing independently checked export correspondence. Do not label a derived estimate tree as original Kraken read classification.

FAPROTAX continues to select Bacteria/Archaea only. For these two fungal abundance tables, an explicitly requested FAPROTAX module should emit its current `skipped` result/reason and no fabricated functional abundance. The optional package/version preflight still applies unless changed in a separate justified patch. For a synthetic mixed/prokaryotic Kraken case, check percentages/count-basis labels and clarify that functional values derive from estimated counts and taxon-based inference.

For the two-sample descriptive cohort, expect Bray/Jaccard distances, eligible one-dimensional PCoA/PCA scores, composition and shared-taxa outputs. Preserve skips for 2D PCoA/PCA plots, NMDS, PERMANOVA, and dispersion inference where prerequisites fail. Do not manufacture replicates, relax existing gates, infer blank-control status from `Knegative`, or subtract that sample from Biomecon. Use neutral placeholder grouping until the user's experimental metadata supplies meaning.

## 8. Provenance, transactions, and independent verification

Add all consumed reports and per-sample params to config defaults/path resolution, sample-ID validation, initial/full inventories, duplicate-physical-file checks, output-root overlap checks, mutation guards, transaction identity, and manifest input collections. Optional consumed logs/versions/taxonomy mappings need the same treatment. Keep validation read-only and cache-only by default; existing persistent diagnostics, credential sanitization, locks, rollback, and ownership rules remain in force.

Write manifest schema v2 revision 5 for the new population/count-basis contract. Add a strict revision-5 validator and retain readers for revisions 1–4. New manifests must validate classifier-specific fields and equalities, evidence availability, report fingerprints, inactive thresholds, and derived export basis. Keep config schema version 1 for additive optional keys unless the implementation introduces a breaking config change. Record the schema decision in release notes.

`tests/verify_release_run.R` currently hardcodes minimap2, one database set, null workflow metadata, and revisions 3/4. Extend it using independently read, fingerprinted producer inputs and known fixture identities. Do not merely permit either classifier or trust the manifest's declared denominators. Independently recompute N/C/U/E/T/D, Bracken row matching, unchanged matrices, generated kreport root/U/percentages, and viewer correspondence; retain AAy/minimap2-specific assertions for those fixtures. New source files must appear in maintained source inventories.

## 9. Cohort helper: second phase after core support

Implement `analysis/utils/merge_wf16s_cohort.R` after single-sample support passes. Keeping it under `analysis/` lets existing source-provenance scope include it; the draft's top-level `utils/` path would need additional provenance rules.

Define repeated `--input <producer-output-root>` arguments, `--output-dir`, `--project-name`, and `--validate-only`. Use each inner output root. Require compatible classifier, NCBI database resources, rank, preprocessing thresholds, active classifier thresholds, Bracken settings, and available workflow/tool evidence. Record unknown effective model length/database checksums rather than claiming byte-identical upstream databases from matching URLs alone. Reject mixed classifiers, incompatible contracts, duplicate/colliding sample IDs, duplicate taxa, source reuse, unsafe output overlap, and conflicting TaxID/lineage identities.

Join by validated full `tax` paths with stable ordering; preserve original counts, fill absent samples with integer zero, recompute `total` safely, and retain exactly one synthetic unclassified row first. Keep assignment/report mappings pointing to original files. Preserve every original complete params file through `params_by_sample`; the representative params file is not provenance for the whole cohort. Emit a hashed merge manifest with each source abundance/params/report mapping and the compatibility result. Consumer validation must inventory and verify this merge evidence if it is used to assert provenance.

Generate complete metadata with neutral `Group=Unspecified` and exact sample IDs, plus a config using the new mappings, a separate ITS taxonomy cache, ITS plot bounds, and a separate output directory. Clearly label the cohort descriptive and experimental metadata unconfirmed. Do not invent treatment/control labels.

Validate all sources before writing. Stage helper outputs beside the destination and atomically publish only after complete success; refuse an existing destination by default. `--validate-only` performs no filesystem/network/cache mutation. Test failure after preparation leaves no partial final bundle and preserves inputs/unrelated destination contents. No FASTQ copying, classifier rerunning, or sequencing-file merging is needed.

## 10. Implementation sequence and acceptance gates

### Concrete code touch points

| File or area | Required work |
| --- | --- |
| `analysis/utils/io.R` | `read_upstream_params`, `extract_upstream_contract`, both assignment passes, `build_context`, classifier-aware reconciliation/BAM discovery |
| New `analysis/utils/classifier_contract.R` | Strict report parsers, TaxID-backed row mapping, normalized producer comparison, population accounting |
| `analysis/utils/config.R`, `config.example.yml` | New input defaults/paths/mapping validation and documented count semantics |
| `analysis/utils/preflight.R` | Inventories, duplicate sources, overlap guards, resolver classifier/evidence arguments and pre-staging checks |
| `analysis/utils/ncbi_taxonomy.py` | Explicit assignment classifier, parser parity, narrow missing-rank handling and evidence/provenance integration |
| `analysis/00_run_pipeline.R` | Source new adapter, initial/full input guards, summaries, transaction inputs, manifest revision/count-basis payload |
| `analysis/01_qc_diagnostics.R` | Read versus estimate fields, classifier-aware fallback/assertions/plots, inactive alignment diagnostics |
| `analysis/02_alpha_diversity.R`, `04_taxa_composition.R`, `08_faprotax.R` | Audit/revise count labels and exported basis; retain established numeric calculations and eligibility gates |
| `analysis/07_kreport_pavian.R`, `analysis/utils/kreport.R` | Resolver call, report inputs, estimate export basis/provenance and applicable output checks |
| Krona/explorer/Sankey renderers and correspondence verifiers | Propagate basis to payloads/captions without changing counted source matrix or official report format |
| `analysis/utils/manifest.R` | Strict revision-5 fields/equalities and backward readers |
| `tests/verify_release_run.R` | Independent classifier/report/accounting checks and preservation of old fixture-specific expectations |
| `tests/testthat`, `tests/test_ncbi_taxonomy.py`, synthetic fixtures | Parser parity, denominator, identity, zero-case and transaction regressions |
| `.github/workflows/ci.yml` | New Kraken integrations and existing cross-platform regression/export gates |
| New `analysis/utils/merge_wf16s_cohort.R` | Second-phase merger and hashed per-source provenance |
| README, CHANGELOG, VERSION, CITATION/release checks | Supported formats, changed semantics, workflow examples and coherent proposed version |

1. **Contracts and input adapter:** implement parameter validation, config maps, strict R/Python schemas, report readers, row identity and N/C/U/E/T/D accounting. Add adversarial fixtures before relaxing the classifier gate. Establish narrow placeholder rules from producer evidence.
2. **Context, QC, and exports:** wire count basis and raw accounting, guard all BAM paths, fix labels/ranges, preserve matrix algorithms, and propagate estimate provenance into all report viewers.
3. **Manifest and verification:** implement revision 5, fingerprint every consumed source, preserve old revisions, and extend the independent release verifier.
4. **Core validation:** run existing minimap2 regressions plus the synthetic Kraken fixtures on Windows and Ubuntu. Obtain passing core gates before building the merger.
5. **Cohort helper and real validation:** validate/run each supplied sample separately, then their descriptive cohort. Inspect TSV/JSON counts and rendered plots/HTML. Publish a validation record with commands, outcomes, warnings/skips, exact candidate SHA if available, input hashes, and limitations.

### Required synthetic cases

| Case | Required outcome |
| --- | --- |
| Existing five-column minimap2, including C+TaxID0 and bamstats | Same counts, metrics, threshold partitions, and source semantics as baseline |
| Six-column Kraken, E=C and E<C | Exact separate accounting; no forced N=T or C=E |
| Same grand totals but a changed species count | Reject Bracken/abundance row mismatch |
| Missing/extra report, wrong sample or reused physical source | Reject before staging |
| Wrong-classifier schema; mixed widths within/across chunks; trailing tabs | Reject with physical line and sample/file context |
| C+0/U+positive; duplicate read IDs across chunks; bad canonical TaxID/length | Reject in both R and Python |
| Valid A hit token; malformed hit token; paired/protein separators | Accept supported DNA syntax; reject unsupported syntax explicitly |
| Changed assignment/report/params between validation and execution | E_INPUT_CHANGED and preserved final outputs/cache |
| Name alias, placeholder rank, and contradictory named ancestry | Accept evidence-backed identity; reject contradictions/ambiguity |
| Bracken species without direct assignment | Resolve from verified mapping evidence or fail with actionable missing-evidence error |
| Kraken reports with assignments absent | Verify report/matrix accounting; label read diagnostics unavailable |
| All unclassified, zero E with positive raw C, and no-read input | Preserve supported zero-classified module gates; no-read input stays unsupported unless separately specified |
| Fungal-only and mixed/prokaryotic abundance | Fungal FAPROTAX skip; valid inference/count-basis behavior for eligible taxa |
| Two-sample and >=3-sample synthetic cohorts | Real pair obeys existing gates; suitable >=3 design exercises 2D ordination; replicated >=2-group fixture exercises inference gates |
| Derived kreport/Krona/explorer/Sankey enabled and disabled | Correct basis, original report unchanged, correspondence and ownership verified |
| Merger incompatibility, unsafe paths, preexisting output, interrupted publication | Fail cleanly and preserve original inputs and unrelated files |

### Commands and execution evidence

Use the repository's locked environment and existing CI setup. Do not install/update unrelated packages or alter `renv.lock` merely to make a gate pass. Commands for implementation validation include:

```powershell
Rscript tests/testthat.R
python -m unittest discover -s tests -p "test_*.py"
node --test tests/test_taxonomy_sankey_client.mjs
python tests/check_release_metadata.py
git -c safe.directory=D:/W/AAy_Amplicon diff --check
```

Add small tracked fixtures under `tests/fixtures/synthetic_kraken2/`, including a dedicated offline cache and complete config. The current audit helper/evidence can remain handoff material under `patch/`; it is not a substitute for automated integration tests. Do not commit the two full real bundles to get CI inputs.

Add CI integrations using the new fixture config and the independent verifier, both on Windows and Ubuntu, alongside existing minimap2/zero-classified/phylogenetic/FAPROTAX checks. Extend cross-OS canonical Sankey checks for the Kraken estimate basis without expecting machine-specific manifests to be byte-identical. Do not hardcode historical test counts.

For each new fixture and real config, run `--validate-only` and compare before/after input/cache/output inventories (including mtime where the contract requires it). A nonexistent output directory must remain nonexistent. Negative validation must leave no staging/output/cache writes. Then run the maintained module set and `--krona --pavian` as appropriate, with deterministic seeds. Check export-enabled and disabled ownership/overwrite behavior, plus a late-failure rollback case.

Real configs and their actual validation commands are implementation deliverables. They must reference the paths in section 2, both new report mappings, the correct output params, separate output/cache locations, and 250–1,000 bp display limits. Test BOTH supplied samples; Biomecon alone is insufficient. Confirm the exact measured numbers above and the FAPROTAX/ordination skips. Keep private run metadata local. Classifier installation and raw-read reprocessing are not required for this postprocessor test.

When VERSION advances to 0.5.0, update CHANGELOG, README, config example, CITATION and release metadata/fixtures as required by `tests/check_release_metadata.py`. Preserve historical release immutability. Local passes, unavailable gates, and cross-platform CI passes must be reported separately. Do not claim a complete suite pass when it times out or a child process crashes.

### Definition of done

Both real bundles pass strict read-only validation and end-to-end execution with their original species matrix retained. Their raw-read and estimate counts reconcile under the separate models, with no invented counts or silent identity mismatch. Minimap2 behavior remains regression-protected. The descriptive cohort runs with the correct statistical skips. Revision-5 provenance and independent export verification pass. New synthetic CI gates cover valid mismatches, malformed inputs, identity conflicts, zero cases, and transactional failures. The delivery record is concrete enough to reproduce every claimed gate; publication remains governed by the execution chat's explicit authorization.

## 11. What this review actually validated

Completed: local baseline/schema/call-site inspection, full scans of 66,847 assignment rows, raw-report population checks, complete Bracken-to-abundance row checks for both samples, source hashes, alias/placeholder evidence, and examination of existing QC, module gates, provenance and release-verifier assumptions. The audit JSON was regenerated from the included script.

Not performed: pipeline implementation, new-input pipeline execution (the current classifier gate rejects it), taxonomy refresh, full R/Python regression suites for proposed changes, upstream reclassification, cohort generation, commit/push/tag, or release publication. Public-source fetch limitations are stated in section 2 rather than treated as successful verification.
