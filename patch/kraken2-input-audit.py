"""Read-only audit of the two local ITS bundles; JSON goes to stdout.

This is handoff evidence, not a production parser or taxonomy resolver.
No sequences are read and no network requests are made.
"""

import csv
import hashlib
import json
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": path.relative_to(ROOT).as_posix(), "sha256": digest.hexdigest()}


def whole(text):
    value = Decimal(text)
    assert value.is_finite() and value >= 0 and value == int(value), text
    return int(value)


def audit(name):
    run = ROOT / name / name
    output = run / "output"
    assignment_path, = (output / "reads_assignments").glob("*.tsv")
    kraken_path, = (output / "kraken2").glob("*.txt")
    bracken_path, = (output / "bracken").glob("*.report")
    abundance_path = output / "abundance_table_species.tsv"
    params_path = output / "params.json"
    params = json.loads(params_path.read_text(encoding="utf-8"))
    launch_params = json.loads((run / "params.json").read_text(encoding="utf-8"))
    statuses = Counter()
    assigned_taxids = Counter()
    lengths = []
    seen = set()
    taxid_lineages = defaultdict(set)
    with assignment_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            assert len(fields) == 6, line_number
            status, read_id, taxid, length, hits, lineage = fields
            assert status in {"C", "U"} and read_id and read_id not in seen
            assert re.fullmatch(r"0|[1-9][0-9]*", taxid)
            assert (status == "U") == (taxid == "0")
            assert re.fullmatch(r"[1-9][0-9]*", length)
            assert hits and all(re.fullmatch(r"(?:0|[1-9][0-9]*|A):[0-9]+", t)
                                for t in hits.split(" "))
            seen.add(read_id)
            statuses[status] += 1
            assigned_taxids[taxid] += 1
            lengths.append(int(length))
            if status == "C":
                assert len(lineage.split("|")) == 7 and all(lineage.split("|"))
                taxid_lineages[taxid].add(lineage)

    with abundance_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sample, = [key for key in rows[0] if key not in {"tax", "total"}]
    abundance = {}
    unclassified = None
    for row in rows:
        assert len(row["tax"].split(";")) == 8
        count = whole(row[sample])
        assert count == whole(row["total"])
        if row["tax"].startswith("Unclassified;"):
            assert unclassified is None
            unclassified = count
        else:
            assert row["tax"] not in abundance
            abundance[row["tax"]] = count

    with kraken_path.open(encoding="utf-8", newline="") as handle:
        report = list(csv.reader(handle, delimiter="\t"))
    assert all(len(row) == 6 for row in report)
    root_row, = [row for row in report if row[3] == "R" and row[4] == "1"]
    u_row, = [row for row in report if row[3] == "U" and row[4] == "0"]
    assert whole(root_row[1]) == statuses["C"]
    assert whole(u_row[1]) == whole(u_row[2]) == statuses["U"] == unclassified
    report_direct = {row[4]: whole(row[2]) for row in report}
    assert len(report_direct) == len(report)
    assert {taxid: count for taxid, count in report_direct.items() if count} == dict(assigned_taxids)

    with bracken_path.open(encoding="utf-8", newline="") as handle:
        bracken = list(csv.DictReader(handle, delimiter="\t"))
    reconstructed = {}
    aliases = []
    placeholder_differences = []
    for row in bracken:
        assert row["taxonomy_lvl"] == "S"
        taxid = row["taxonomy_id"]
        lineage, = taxid_lineages[taxid]
        parts = lineage.split("|")
        # Explicit local comparison only: observed 8-rank tables add Fungi.
        def rank_matches(left, right, rank):
            missing = re.fullmatch(r"unclassified (.+) " + rank, right)
            return left == right or (missing is not None
                                     and left == missing[1] + "_Incertae_sedis")

        candidates = [path for path in abundance
                      if all(rank_matches(left, right, rank) for left, right, rank in zip(
                          path.split(";")[:1] + path.split(";")[2:], parts,
                          ("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
                          strict=True))]
        path, = candidates
        if path.split(";")[:1] + path.split(";")[2:] != parts:
            placeholder_differences.append({"taxid": taxid, "abundance_path": path,
                                           "assignment_lineage": lineage})
        assert path not in reconstructed
        count = whole(row["new_est_reads"])
        assert whole(row["kraken_assigned_reads"]) + whole(row["added_reads"]) == count
        reconstructed[path] = count
        if row["name"] != parts[-1]:
            aliases.append({"taxid": taxid, "bracken_name": row["name"],
                            "lineage_name": parts[-1], "estimated_count": count})
    assert reconstructed == abundance
    estimates = sum(abundance.values())
    assert estimates <= statuses["C"]
    stdout = (run / "nextflow.stdout").read_text(encoding="utf-8")
    versions = sorted(set(re.findall(r"wf-16s v([0-9.]+)", stdout)))
    revisions = sorted(set(re.findall(r"revision: ([0-9a-f]+)", stdout)))
    return {
        "sample": sample, "workflow_versions_from_log": versions,
        "workflow_revisions_from_log": revisions,
        "params": {key: params.get(key) for key in (
            "classifier", "database_set", "taxonomic_rank", "min_len", "max_len",
            "min_read_qual", "min_percent_identity", "min_ref_coverage",
            "abundance_threshold", "kraken2_confidence", "bracken_threshold", "bracken_length")},
        "launch_params_missing_contract_fields": sorted(
            key for key in ("database_sets", "bracken_length") if key not in launch_params),
        "assignment_rows": len(seen), "raw_classified": statuses["C"],
        "raw_unclassified": statuses["U"], "c_taxid_zero": 0,
        "length_range": [min(lengths), max(lengths)],
        "abundance_rows": len(rows), "species_estimate_sum": estimates,
        "abundance_unclassified": unclassified,
        "abundance_table_sum": estimates + unclassified,
        "raw_classified_minus_species_estimates": statuses["C"] - estimates,
        "positive_abundance_domains": sorted({path.split(";")[0]
                                              for path, count in abundance.items() if count}),
        "bracken_taxid_backed_row_match": True, "name_aliases": aliases,
        "kraken_direct_counts_match_assignments": True,
        "placeholder_differences": placeholder_differences,
        "files": [fingerprint(path) for path in (
            params_path, abundance_path, assignment_path, kraken_path, bracken_path,
            output / "versions_all.txt", run / "nextflow.stdout")],
    }


if __name__ == "__main__":
    result = {"audit_date": "2026-10-10", "method": "local files, no sequence or network access",
              "samples": [audit(f"wf-16s_{name}-ITS-kraken2")
                          for name in ("Biomecon", "Knegative")]}
    print(json.dumps(result, ensure_ascii=False, indent=2))
