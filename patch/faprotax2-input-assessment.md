# FAPROTAX2 local sequence-input assessment

**Scope:** Read-only inventory of `other-references-inputs` for Phase 1. Original inputs were not modified, copied, or committed. This report contains aggregate metadata only; an ignored local inventory retains exact paths and per-file SHA-256 hashes for reproducibility.

## Findings

The inventory contained 1,438 files and 12 FASTA files. Six FASTAs named `expr.asv.fasta` had matching `expr.otu_table.txt` feature tables beside their analysis outputs. All six are 16S: three V3–V4 datasets and three archaeal-targeted V4–V5 datasets. The associated project metadata names Illumina and DADA2, so these are Illumina-derived representative ASVs rather than ONT reads. No matched ITS dataset was found.

For every pair, the FASTA first-token identifiers and table feature identifiers matched exactly and in the same order. Neither side had duplicate feature identifiers. The corresponding `expr.tax_table.txt` files also covered the same feature IDs exactly, with seven taxonomy columns. The count tables store samples in rows and features in columns. Every value parsed as finite, nonnegative, and integer-valued; every feature had a positive total across the listed samples. No sample labels, feature identifiers, taxonomy strings, or sequence bases are reproduced here.

| Local candidate | Marker | Features | Samples | FASTA length, min / median / max (nt) | Exact FASTA-to-table IDs |
|---|---|---:|---:|---:|---|
| A | 16S V3–V4 | 686 | 4 | 204 / 425 / 430 | Yes, same order |
| B | 16S V4–V5, archaeal-targeted | 66 | 4 | 385 / 385 / 386 | Yes, same order |
| C | 16S V3–V4 | 1,315 | 7 | 202 / 424 / 431 | Yes, same order |
| D | 16S V4–V5, archaeal-targeted | 357 | 7 | 230 / 385 / 386 | Yes, same order |
| E | 16S V3–V4 | 602 | 3 | 334 / 428 / 431 | Yes, same order |
| F | 16S V4–V5, archaeal-targeted | 44 | 3 | 385 / 385 / 387 | Yes, same order |

All six FASTAs had unique nonempty first-token IDs, nonempty sequences, and no ambiguous bases or gap/invalid symbols under a case-insensitive IUPAC check. Sequence lengths range from 202 to 431 nt. The smallest complete dataset is candidate F (44 ASVs and 3 samples); candidate A is a larger bacterial V3–V4 option (686 ASVs and 4 samples). These are whole existing feature sets, not sequences synthesized from taxonomy labels or subsets reconstructed from species abundance tables.

## Use and limits

Candidate F is the smallest real 16S input pair for sequence-mode testing if the official mapper supports this archaeal-targeted V4–V5 marker. Candidate A is available if the mapper requires a bacterial V3–V4 set. The source metadata confirms the Illumina/DADA2 context; it does not establish that either marker span is supported by FAPROTAX2. Use an official mapper contract and a successful native-Windows run to settle that question before calling either pair a valid end-to-end sequence-mode fixture.

The table values are integer-valued, but this inventory does not establish their producer's precise abundance semantics or compatibility with a particular mapper input format. The files were left in their existing samples-by-features orientation. No transpose, filtering, ID rewriting, or count normalization was performed. The local ignored evidence is at `output/faprotax2-audit/inputs/inventory.private.json`; it must remain local because it contains absolute input paths and hashes.
