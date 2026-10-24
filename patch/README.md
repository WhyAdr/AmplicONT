# Plans and validation evidence

Implementation plans and review records live here so the repository root stays
focused on pipeline configuration and project metadata. Historical plans describe
their original scope; the current implementation and usage are documented in the
[project README](../README.md).

## FAPROTAX2 migration

- [Native Windows assessment and legacy baseline](faprotax2-migration-assessment.md)
- [Official mapper contract and provenance](faprotax2-upstream-contract.md)
- [Architecture integration audit](faprotax2-architecture-assessment.md)
- [Representative 16S input assessment](faprotax2-input-assessment.md)
- [Taxonomy proof of concept](faprotax2-poc-results.md)
- [Production taxonomy integration and native validation](faprotax2-production-validation.md)
- [Sequence feasibility and native portability trial](faprotax2-sequence-feasibility.md)
- [Final native sequence trial and blocked gate](faprotax2-native-sequence-trial.md)
- [Independent two-backend comparison and remaining limits](faprotax2-validation-report.md)

## Version 0.5.0: Kraken2 and Bracken

- [Original integration draft](kraken2-support-implement.md)
- [Refined implementation handoff](kraken2-support-implementation-refined.md)
- [Read-only input audit](kraken2-input-audit.py) and [audit evidence](kraken2-input-audit.json)
- [Implementation review and local validation](v0.5.0-review-validation.md)
  and [machine-readable validation evidence](v0.5.0-review-validation.json)

Run the audit from the repository root with `python patch/kraken2-input-audit.py`.
It requires the original local ITS bundles, which are ignored by Git.

## Version 0.4.9

- [NCBI taxonomy error and rank plan](v0.4.9-ncbi-taxonomy-error-rank-patch-plan.md)
- [Validated patch plan](ont-wf16s-v0.4.9-astra-validated-patch.md)
- [Release hardening follow-up](releasev0.4.9-hardening-2.md)
- [Release validation checklist](v0.4.9-release-validation-checklist.md)

## Earlier versions and integrations

- [Version 0.4.8 Sankey renderer plan](wf16s-postprocess-v0.4.8-sankey-renderer-plan-refined.md)
- [Version 0.4.7 refined patch plan](wf16s-postprocess-v0.4.7-patch-refined.md)
- [Version 0.4.6 alpha diversity and audit plan](v0.4.6-alpha-diversity-and-audit-implementation-plan.md)
- [FAPROTAX integration plan](faprotax1.2.12-integration-plan.md)
  and [concordance review](faprotax1.2.12-concordance-review.md)
- [Krona integration plan](krona-script-integration-plan.md)

Older patch plans and release audits remain alongside these indexed documents.
