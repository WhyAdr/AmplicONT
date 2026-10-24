# Native Windows FAPROTAX2 sequence trial

**Status: blocked; sequence mode is not verified.** The final one-feature native Windows run completed MAFFT alignment and generated an APPLES jplace file, then exited with status 1 at FAPROTAX's reference-tree consistency check. It produced no functional abundance table or annotation report. Alignment and jplace generation do not count as sequence-mode success. No production adapter was integrated; P4 remains unimplemented.

This report records a single bounded trial using one sample-derived 16S feature, its exactly matching representative sequence, and three abundance-table sample columns. Feature IDs and sample labels are omitted. Input sequences, raw tables, and private labels remain local and ignored.

## Trial result

The isolated run used the complete official FAPROTAX2 reference package copied to a disposable run directory. The original compressed reference alignment had 26,632 records while the official reference tree had 13,816 tips. A streaming projection retained exactly one record for every tree tip, rejected duplicate or missing tip IDs, excluded 12,816 alignment records absent from the tree, and verified that the retained sequence bases were unchanged.

MAFFT 7.526 completed the fragment alignment. APPLES 2.0.11 then completed its default FastTree re-estimation and wrote a jplace file for the one query. The mapper stopped before function mapping with:

> ERROR: The reference tree in the placement file differs from the FAPROTAX reference tree.

The mapper exit status was 1. No functional table or report was written. The jplace file contains one placement query; it is evidence that placement ran, not evidence that FAPROTAX sequence mode succeeded.

### Direct reference-tree comparison

The jplace Newick tree was compared with the original tree using exact tip-set and descendant-tip bitset signatures after removing jplace edge-number annotations. The disposable tree was byte-identical to the official reference tree before APPLES ran.

| Check | Result |
|---|---:|
| Unique tip labels in each tree | 13,816; identical sets |
| Root child counts | 2 and 2 |
| Root bipartition | Identical |
| Internal unrooted splits | 13,813 per tree |
| Shared internal splits | 13,798 |
| Reference-only / jplace-only internal splits | 15 / 15 |
| Shared splits whose branch-length sums differ | 27,614 / 27,614 |
| Median / maximum absolute branch-length-sum difference | 0.029209 / 3.425421 |

The reference and jplace trees therefore differ in 15 internal resolutions as well as branch lengths. The roots and tip labels match, so a label or root-format adjustment cannot repair this check. Bypassing the tree comparison or forcing placements onto the original tree would change placement semantics and was not done. APPLES ran with one thread and without `-D`, preserving its default re-estimation behavior for the ML-derived reference tree.

## Native runtime and provenance

| Component | Version/artifact | SHA-256 | License evidence |
|---|---|---|---|
| MAFFT native Windows package | `mafft-7.526-win64-signed.zip`, 27,476,513 bytes | `7C3A508C8512AF04007D68019C6FBD4FF00D891244B58F87EC9AA0E53FBAE8A7` | Upstream identifies MAFFT as GPL; see the [official Windows package page](https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html). |
| APPLES source distribution | `apples-2.0.11.tar.gz`, 760,816 bytes | `0df2367ceac167c6c347954911643a0305c2f293a5f31a76bd86bd568bf12115` | The PyPI sdist metadata says `License: UNKNOWN`; the [official v2.0.11 source tag](https://github.com/balabanmetin/apples/tree/v2.0.11) identifies GPL-3.0. |

The isolated runtime used CPython 3.12.10 x64 in `output/faprotax2-runtime/sequence-py312-clean/` and R 4.5.3 with an ignored R user library. The APPLES wheel dependencies were downloaded at pinned versions and their archive hashes matched PyPI metadata. The isolated Python package inventory was:

```text
apples==2.0.11
niemads==1.0.16
numpy==1.26.4
setuptools==75.8.0
TreeCluster==1.0.5
TreeSwift==1.1.50
wheel==0.45.1
```

`pip check` reported no broken requirements. APPLES' bundled FastTree reported version 2.1.11. The downloaded and unpacked upstream archives were kept unchanged; all mapper and APPLES compatibility edits were applied to scratch copies only.

The official reference tree SHA-256 is `1b40b2e3f4165b6eea49290064aadaff3c74b4744c16e91c35506afcfdfa977d`; the disposable v4 tree had the same SHA-256. Private before/after manifests for the official unpacked source and reference files show no changes to those originals.

## Windows fixes exercised in scratch copies

These edits demonstrate which native portability barriers can be passed; none has been integrated into production:

1. Corrected the Windows decompression output target so the decompressed file is written to the intended uncompressed output path rather than overwriting the `.gz` input.
2. Invoked MAFFT's Windows batch launcher through `cmd.exe` with explicit arguments, placed `--out` before `--addfragments`, captured its status, and rejected a missing or empty output.
3. Invoked APPLES with the isolated venv's Python executable and an argument vector, then checked its exit status and jplace output.
4. Guarded APPLES' POSIX-only `fork` start-method call. Windows keeps APPLES' serial query-processing branch.
5. Changed APPLES' Windows TreeCluster subprocess to use the venv interpreter and script path, with explicit exit-status handling.
6. Kept the representative-sequence preprocessing path serial on Windows, preserving input order and leaving the POSIX pool path intact.
7. Decompressed and streamed the official alignment through an exact tree-tip projection before APPLES. This preserved every retained sequence and the reference tree, while excluding alignment records that had no tree tip.

These adaptations let the one-query run reach jplace creation. They did not make the changed APPLES backbone compatible with the FAPROTAX reference-tree gate.

## Resource limits and remaining work

The host exposed 8,277,602,304 bytes of RAM (about 7.71 GiB). During the one-thread full-reference FastTree step, observed process working sets were about 1.13–1.19 GiB. Available host memory fell to about 444 MiB in the final run and about 116 MiB in an earlier scratch attempt. Only one feature was run; no 44-feature run was launched. This is insufficient evidence for throughput or multi-feature memory requirements.

The sequence-mode phase remains blocked until an authorized, scientifically valid route produces a jplace tree compatible with the FAPROTAX reference and completes functional annotation end to end. No tree comparison was bypassed, no reference topology or branch lengths were overwritten, and no alternate success was inferred from the partial outputs. Taxonomy-mode validation is separate.

## Local reproducibility artifacts

The following ignored, checkout-relative artifacts preserve the final scratch implementation and aggregate evidence for local review:

- `output/faprotax2-audit/sequence/mapper-native-windows-v4.patch`
- `output/faprotax2-audit/sequence/patch_mapper_native_windows_v4.py`
- `output/faprotax2-audit/sequence/mapper-scratch-v4/faprotax.R`
- `output/faprotax2-runtime/apples-native-windows-v4.patch`
- `output/faprotax2-runtime/patch_apples_native_windows_v4.py`
- `output/faprotax2-runtime/project_reference_to_tree_tips.py`
- `output/faprotax2-audit/sequence/run_native_single_matched_v4.ps1`
- `output/faprotax2-audit/sequence/prepare_single_matched_v4.py`
- `output/faprotax2-audit/sequence/native-patched-single-matched-v4/tree-comparison.private.json`
- `output/faprotax2-runtime/pypi-distribution-verification.private.json`

The runner and fixture-preparation script are local-only and resolve private inputs from this checkout; their machine-specific values are not reproduced here. Private FASTA/table fixtures, the jplace file, and detailed run logs stay under ignored `output/faprotax2-audit/sequence/`; they are not included here or committed. The ignored artifacts are available only in the current local checkout.
