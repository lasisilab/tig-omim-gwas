#!/usr/bin/env python3
"""Reproduce the OMIM analysis from literature tables and saved API responses.

Run with Python 3.11+ from any directory. This script makes no network requests,
does not read previously derived CSVs, and leaves inputs unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "scripts"))

import pandas as pd

import omim_ids as ids
import omim_gene_maps as gene_maps


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def csv_bytes(frame):
    # The archived responses were originally processed with pandas' LF CSV output.
    return frame.to_csv(index=False, lineterminator="\n").encode("utf-8")


def csv_sha256(frame):
    return hashlib.sha256(csv_bytes(frame)).hexdigest()


def contained_path(root, relative):
    path = (root / relative).resolve()
    require(path.is_relative_to(root.resolve()), "Input path is outside its expected folder")
    return path


def verify_inputs(inputs):
    records = read_json(inputs / "checksums.json")
    require(len(records) == len({row["file"] for row in records}), "Duplicate input checksum paths")
    verified = []
    for record in records:
        path = contained_path(inputs, record["file"])
        actual = ids.sha256(path)
        require(actual == record["sha256"], f"Input checksum changed: {record['file']}")
        verified.append({"file": record["file"], "sha256": actual,
                         "bytes": path.stat().st_size, "verified": True})
    inventory = {p.relative_to(inputs).as_posix() for p in inputs.rglob("*")
                 if p.is_file() and p.name != "checksums.json"}
    require(inventory == {row["file"] for row in records}, "Input checksum inventory is incomplete")
    return pd.DataFrame(verified)


def reconstruct_names(inputs, unique):
    pointer = read_json(inputs / "local-data/omim/latest.json")
    directory = contained_path(inputs, pointer["run_directory"])
    require(directory.is_relative_to(inputs / "local-data/omim"), "Unexpected names run directory")
    run = read_json(directory / "run.json")
    require(csv_sha256(unique) == run["seed_sha256"], "Recomputed cited ID checksum changed")
    requests = read_json(directory / "requests.json")
    entries = {}
    for request in requests:
        response = contained_path(directory, request["response_file"])
        require(ids.sha256(response) == request["sha256"], "Name response checksum changed")
        requested = set(request["requested_ids"])
        returned = []
        for wrapper in read_json(response)["omim"]["entryList"]:
            entry = wrapper["entry"]
            identifier = str(entry["mimNumber"])
            require(identifier in requested and identifier not in entries,
                    "Unexpected or duplicate ID in a name response")
            prefix = entry.get("prefix", "")
            title = entry.get("titles", {}).get("preferredTitle", "")
            require(prefix in ids.ENTRY_TYPES and bool(title), "Unrecognized OMIM prefix or missing title")
            entries[identifier] = {
                "omim_id": identifier, "preferred_title": title, "prefix": prefix,
                "entry_type": ids.ENTRY_TYPES[prefix], "status": entry.get("status", ""),
                "alternative_titles": entry.get("titles", {}).get("alternativeTitles", ""),
                "included_titles": entry.get("titles", {}).get("includedTitles", ""),
                "retrieved_at_utc": request["retrieved_at_utc"],
                "omim_url": f"https://omim.org/entry/{identifier}",
            }
            returned.append(identifier)
        require(len(returned) == len(request["requested_ids"]) and set(returned) == requested,
                "Name response does not cover its requested IDs")
        require(returned == request["returned_ids"], "Name response differs from returned-ID log")
    require(set(entries) == set(unique["omim_id"]), "Name responses do not cover all cited IDs")
    names = pd.DataFrame([entries[identifier] for identifier in unique["omim_id"]])
    require(len(names) == run["unique_ids"], "Name metadata count changed")
    require(csv_sha256(names) == run["names_sha256"], "Names rebuilt from responses differ from archived checksum")
    return names, ids.resolve_moved_titles(names), directory, requests


def reconstruct_maps(inputs, unique, names):
    pointer = read_json(inputs / "local-data/omim-genes/latest.json")
    directory = contained_path(inputs, pointer["run"])
    require(directory.is_relative_to(inputs / "local-data/omim-genes"), "Unexpected mappings run directory")
    require(ids.sha256(directory / "run.json") == pointer["run_sha256"], "Mapping run checksum changed")
    run = read_json(directory / "run.json")
    require(csv_sha256(unique) == run["seed_sha256"], "Mapping seed checksum changed")
    require(sorted(set(names["resolved_omim_id"])) == run["requested_ids"], "Mapping targets changed")
    require(ids.sha256(inputs / "data/omim/gene-query.json") == run["query_sha256"],
            "Public gene query differs from the archived query")
    _, mappings, coverage, requests = gene_maps.read_maps(directory)
    require(set(coverage["resolved_omim_id"]) == set(names["resolved_omim_id"]), "Mapping coverage mismatch")
    return mappings, coverage, directory, requests


def make_pipeline(cells, mentions, names, coverage, evidence):
    rows = []
    for group in ("all", "pigmentation", "hair"):
        source = cells if group == "all" else cells.loc[cells["group"].eq(group)]
        extracted = mentions if group == "all" else mentions.loc[mentions["group"].eq(group)]
        identifiers = set(extracted["omim_id"])
        returned = names.loc[names["omim_id"].isin(identifiers)]
        targets = set(returned["resolved_omim_id"])
        covered = coverage.loc[coverage["resolved_omim_id"].isin(targets)]
        linked = evidence if group == "all" else evidence.loc[evidence["group"].eq(group)]
        with_ids = len(extracted[["source_file", "data_row"]].drop_duplicates())
        moved = int(returned["omim_id"].ne(returned["resolved_omim_id"]).sum())
        symbol_targets = set(covered.loc[covered["approved_genes"].gt(0), "resolved_omim_id"])
        stages = [
            ("source_cells", "Source OMIM cells", "table cells", len(source), "Original selected OMIM columns"),
            ("cells_with_ids", "Cells containing IDs", "table cells", with_ids,
             f"{len(source)-with_ids:,} cells without IDs omitted"),
            ("id_mentions", "Extracted ID mentions", "ID mentions", len(extracted), "Split cells containing multiple IDs"),
            ("cited_ids", "Deduplicate cited IDs", "OMIM IDs", len(identifiers),
             f"{len(extracted)-len(identifiers):,} repeated mentions collapsed"),
            ("returned_names", "Retrieve entry names", "OMIM IDs", returned["omim_id"].nunique(),
             f"{len(identifiers)-returned['omim_id'].nunique():,} cited IDs omitted by API"),
            ("resolved_targets", "Resolve moved IDs", "lookup IDs", len(targets),
             f"{moved:,} moved {'ID' if moved == 1 else 'IDs'}; original citations preserved"),
            ("mapped_targets", "Read approved-symbol mappings", "lookup IDs", len(symbol_targets),
             f"{len(targets-symbol_targets):,} IDs without approved symbols"),
            ("gene_evidence", "Join mappings to cited IDs", "gene–ID rows", len(linked),
             "Keep source ID and trait-group provenance"),
            ("genes", "Deduplicate approved symbols", "gene names", linked["gene_symbol"].nunique(),
             "Exact symbols within each trait group"),
        ]
        for stage, label, unit, count, detail in stages:
            rows.append(dict(group=group, stage=stage, label=label, unit=unit,
                             count=int(count), detail=detail))
    return pd.DataFrame(rows)


def request_provenance(inputs, source, directory, requests):
    return [{"lookup": source, "batch": number, "url": request["url"],
             "requested_ids": "; ".join(request["requested_ids"]),
             "returned_ids": "; ".join(request["returned_ids"]),
             "retrieved_at_utc": request["retrieved_at_utc"],
             "response_file": (directory / request["response_file"]).relative_to(inputs).as_posix(),
             "sha256": request["sha256"]}
            for number, request in enumerate(requests, 1)]


def reproduce(output, inputs=ROOT / "inputs"):
    inputs = inputs.resolve()
    output = output.resolve()
    require(not output.is_relative_to(inputs), "Outputs cannot be written inside inputs")
    verified = verify_inputs(inputs)
    cells = ids.source_id_cells(inputs)
    mentions = ids.extract_ids(cells)
    unique = ids.deduplicate_ids(mentions)
    raw_names, names, names_directory, name_requests = reconstruct_names(inputs, unique)
    mappings, coverage, maps_directory, mapping_requests = reconstruct_maps(inputs, unique, names)
    evidence, all_ids = gene_maps.link_sources(mentions, names, mappings, coverage)
    genes = gene_maps.summarize_genes(evidence)
    pipeline = make_pipeline(cells, mentions, names, coverage, evidence)
    counts = {group: {row.stage: int(row.count) for row in frame.itertuples()}
              for group, frame in pipeline.groupby("group")}

    # Expected values are consulted only after all tables and counts are computed.
    reference_present = (inputs / "expected-counts.json").exists()
    expected = read_json(inputs / "expected-counts.json") if reference_present else {}
    differences = [{"group": group, "stage": stage, "expected": count,
                    "observed": counts.get(group, {}).get(stage)}
                   for group, stages in expected.items() for stage, count in stages.items()
                   if counts.get(group, {}).get(stage) != count]
    if reference_present:
        require(set(expected) == set(counts) and all(set(expected[g]) == set(counts[g]) for g in counts),
                "Expected count dimensions differ from computed stages")
    frames = {
        "omim_source_cells.csv": cells, "omim_id_mentions.csv": mentions, "omim_ids.csv": unique,
        "entry-names.csv": raw_names, "resolved-entry-names.csv": names,
        "gene-mappings.csv": mappings, "mapping-coverage.csv": coverage,
        "id-coverage.csv": all_ids, "gene-evidence.csv": evidence, "genes.csv": genes,
        "omim-retention-pipeline.csv": pipeline, "verified-inputs.csv": verified,
        "api-requests.csv": pd.DataFrame(
            request_provenance(inputs, "entry names", names_directory, name_requests) +
            request_provenance(inputs, "gene maps", maps_directory, mapping_requests)),
    }
    for group in ("pigmentation", "hair"):
        frames[f"{group}-conditions.csv"] = ids.group_table(mentions, names, group)
        frames[f"{group}-genes.csv"] = genes.loc[genes["group"].eq(group)]
    stages = [
        (1, "source_cells", "data/darcy2023_tableS1.xlsx; data/needle2025_table2.csv",
         "Select cited OMIM columns and classify trait groups", "omim_source_cells.csv"),
        (2, "cells_with_ids", "omim_source_cells.csv", "Find cells containing six-digit OMIM IDs", "omim_id_mentions.csv"),
        (3, "id_mentions", "omim_source_cells.csv", "Extract each ID and retain source provenance", "omim_id_mentions.csv"),
        (4, "cited_ids", "omim_id_mentions.csv", "Deduplicate cited OMIM IDs", "omim_ids.csv"),
        (5, "returned_names", "omim_ids.csv; saved names response JSONs", "Rebuild names and verify archived CSV checksum", "entry-names.csv"),
        (6, "resolved_targets", "entry-names.csv", "Resolve moved entries and preserve cited IDs", "resolved-entry-names.csv"),
        (7, "mapped_targets", "resolved-entry-names.csv; saved mapping response JSONs", "Extract approved symbols and verify coverage", "gene-mappings.csv; mapping-coverage.csv"),
        (8, "gene_evidence", "omim_id_mentions.csv; resolved-entry-names.csv; gene-mappings.csv", "Join mappings to source IDs and groups", "gene-evidence.csv; id-coverage.csv"),
        (9, "genes", "gene-evidence.csv", "Deduplicate exact approved symbols per trait group", "genes.csv; pigmentation-genes.csv; hair-genes.csv"),
    ]
    stage_frame = pd.DataFrame(stages, columns=["step", "stage", "inputs", "operation", "outputs"])
    frames["stage-input-output.csv"] = pipeline.merge(stage_frame, on="stage", validate="many_to_one")
    output.mkdir(parents=True, exist_ok=True)
    for filename, frame in frames.items():
        (output / filename).write_bytes(csv_bytes(frame))
    ids.write_json(output / "figure-counts.json", counts)
    software = {"python": platform.python_version(), "pandas": pd.__version__,
                "openpyxl": importlib.metadata.version("openpyxl"), "platform": platform.platform()}
    ids.write_json(output / "software-versions.json", software)
    verification = {
        "success": not differences, "mode": "offline; literature tables and archived raw API responses",
        "input_files_verified": len(verified), "input_checksums_sha256": ids.sha256(inputs / "checksums.json"),
        "names_rebuilt_from_raw_responses": True, "archived_names_checksum_matched": True,
        "reference_present": reference_present,
        "expected_counts_match": not differences if reference_present else None, "differences": differences,
        "cited_ids": len(unique), "resolved_targets": len(coverage), "mapping_rows": len(mappings),
        "gene_evidence_rows": len(evidence), "unique_genes": int(evidence["gene_symbol"].nunique()),
        "genes_by_group": {group: int(count) for group, count in genes.groupby("group").size().items()},
        "source_runs": {"names": names_directory.relative_to(inputs).as_posix(),
                        "gene_maps": maps_directory.relative_to(inputs).as_posix()},
        "output_sha256": {name: ids.sha256(output / name) for name in sorted(frames)},
        "software": software,
    }
    ids.write_json(output / "verification.json", verification)
    print(pipeline[["group", "stage", "count", "unit"]].to_string(index=False))
    print(f"\nVerified {len(verified)} input files; rebuilt names directly from saved responses.")
    print(f"Outputs: {output}")
    require(not differences, "Computed stage counts differ; inspect outputs/verification.json")
    print(f"All {sum(len(stages) for stages in expected.values())} expected stage counts match." if reference_present else
          "All input and response validations passed; no expected-count reference was supplied.")
    return verification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=ROOT / "inputs", help="Input folder (default: package inputs/)")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs", help="Output folder (default: package outputs/)")
    arguments = parser.parse_args()
    reproduce(arguments.output, arguments.inputs)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, AssertionError, json.JSONDecodeError) as error:
        raise SystemExit(f"Reproduction failed: {error}") from None
