# FAPROTAX2 taxonomy POC fixtures

These tables are synthetic and contain no upstream mapper/reference assets or user data. The three full taxonomy paths in `taxonomy_known_species.tsv` correspond to named taxa present in the official mapper reference. Counts and sample labels are test values.

The tables follow the official classical TSV contract: feature IDs in `OTU`, integer abundance columns by sample, and a taxonomy column selected with `-t column:taxonomy`. `taxonomy_unknown_root.tsv` is intentionally outside the reference taxonomy and probes how the upstream mapper handles an unmatched lineage; it is not asserted to be rejected or called unmapped by the adapter.
