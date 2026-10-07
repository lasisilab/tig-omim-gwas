"""Retrieve structured gene mappings for the cited OMIM IDs.

python scripts/omim_gene_maps.py --fetch
Responses and derived OMIM tables stay in ignored local-data/omim-genes/.
Rendering and runs without --fetch only read saved responses.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import time

import pandas as pd

from omim_ids import API, compile_ids, load_names, request_entries, sha256, write_json

ROOT = Path(__file__).resolve().parents[1]
MAP_COLUMNS = ["resolved_omim_id", "condition_name", "entry_prefix", "mapping_index",
    "gene_symbol", "gene_mim_id", "gene_symbols_raw", "approved_symbols_raw",
    "ncbi_gene_ids_raw", "ensembl_ids_raw", "phenotype_mim_id", "phenotype_name",
    "mapping_key", "inheritance", "mapping_field", "mapping_status", "response_file"]


def seed_targets(root=ROOT):
    names, _ = load_names(root)
    if names is None:
        raise ValueError("Retrieve OMIM entry names before gene mappings")
    ids = pd.read_csv(root / "data/omim/omim_ids.csv", dtype=str, keep_default_na=False)
    assert set(ids["omim_id"]) == set(names["omim_id"])
    return names, sorted(set(names["resolved_omim_id"]))


def fetch_maps(root=ROOT):
    config_path = root / "data/omim/gene-query.json"
    config = json.loads(config_path.read_text())
    assert config["api"] == API and config["include"] == "geneMap"
    assert 1 <= config["batch_size"] <= 20
    _, targets = seed_targets(root)
    key = os.environ.get("OMIM_API_KEY", "").strip()
    if not key or key.startswith("op://"):
        raise ValueError("Resolve OMIM_API_KEY before --fetch")
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    directory = root / "local-data/omim-genes" / run_id
    directory.mkdir(parents=True, exist_ok=False)
    (directory / "query.json").write_bytes(config_path.read_bytes())
    requests, seen = [], set()
    for offset in range(0, len(targets), config["batch_size"]):
        batch = targets[offset:offset + config["batch_size"]]
        payload, url = request_entries(batch, key, include=config["include"])
        file = f"response-{len(requests) + 1:02d}.json"
        write_json(directory / file, payload)
        returned = [str(item["entry"]["mimNumber"]) for item in payload["omim"]["entryList"]]
        assert set(returned) == set(batch) and len(returned) == len(batch)
        assert not (seen & set(returned))
        seen.update(returned)
        requests.append({"url": url, "requested_ids": batch, "returned_ids": returned,
            "response_file": file, "sha256": sha256(directory / file),
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat()})
        write_json(directory / "requests.json", requests)
        print(f"Retrieved gene mappings: {len(seen)}/{len(targets)} OMIM entries", flush=True)
        time.sleep(0.7)
    run = {"api": API, "include": config["include"], "requested_ids": targets,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "seed_sha256": sha256(root / config["seed_file"]),
        "query_sha256": sha256(directory / "query.json"),
        "requests_sha256": sha256(directory / "requests.json")}
    write_json(directory / "run.json", run)
    read_maps(directory)
    write_json(root / "local-data/omim-genes/latest.json", {
        "run": directory.relative_to(root).as_posix(),
        "run_sha256": sha256(directory / "run.json")})
    return directory


def flatten_entries(entries):
    """One row per returned mapping and approved symbol; never treat aliases as genes."""
    rows, coverage = [], []
    for entry in entries:
        identifier = str(entry["mimNumber"])
        wrappers = entry.get("phenotypeMapList") or []
        mappings = [(item["phenotypeMap"], "entry.phenotypeMapList") for item in wrappers]
        # A directly queried gene entry is a gene identity, not all its phenotypes.
        if not mappings and entry.get("prefix") in {"*", "+"} and entry.get("geneMap"):
            mappings = [(entry["geneMap"], "entry.geneMap")]
        symbols = set()
        for index, (mapping, field) in enumerate(mappings):
            if field == "entry.phenotypeMapList":
                assert str(mapping["phenotypeMimNumber"]) == identifier
            approved = mapping.get("approvedGeneSymbols") or mapping.get("approvedGeneSymbol") or ""
            assert isinstance(approved, str)
            gene_symbols = sorted(set(s.strip() for s in re.split(r"[,;]", approved) if s.strip()))
            symbols.update(gene_symbols)
            for symbol in gene_symbols or [""]:
                rows.append({"resolved_omim_id": identifier,
                    "condition_name": entry["titles"]["preferredTitle"],
                    "entry_prefix": entry.get("prefix", ""), "mapping_index": index,
                    "gene_symbol": symbol, "gene_mim_id": str(mapping.get("mimNumber") or ""),
                    "gene_symbols_raw": mapping.get("geneSymbols") or "",
                    "approved_symbols_raw": approved,
                    "ncbi_gene_ids_raw": str(mapping.get("geneIDs") or mapping.get("entrezGeneID") or ""),
                    "ensembl_ids_raw": mapping.get("ensemblIDs") or mapping.get("ensemblGeneID") or "",
                    "phenotype_mim_id": str(mapping.get("phenotypeMimNumber") or ""),
                    "phenotype_name": mapping.get("phenotype") or "",
                    "mapping_key": str(mapping.get("phenotypeMappingKey") or ""),
                    "inheritance": mapping.get("phenotypeInheritance") or "",
                    "mapping_field": field,
                    "mapping_status": "approved symbol returned" if symbol else "no approved symbol returned",
                    "response_file": entry["response_file"]})
        coverage.append({"resolved_omim_id": identifier,
            "condition_name": entry["titles"]["preferredTitle"],
            "entry_prefix": entry.get("prefix", ""), "mapping_rows": len(mappings),
            "approved_genes": len(symbols), "gene_symbols": "; ".join(sorted(symbols)),
            "mapping_status": "approved symbol returned" if symbols else
                "no approved symbol returned" if mappings else "no gene mapping returned",
            "response_file": entry["response_file"]})
    maps = pd.DataFrame(rows, columns=MAP_COLUMNS)
    assert not maps.duplicated(["resolved_omim_id", "mapping_index", "gene_symbol"]).any()
    return maps, pd.DataFrame(coverage).sort_values("resolved_omim_id").reset_index(drop=True)


def read_maps(directory):
    run = json.loads((directory / "run.json").read_text())
    assert sha256(directory / "query.json") == run["query_sha256"]
    assert sha256(directory / "requests.json") == run["requests_sha256"]
    entries, requests = [], json.loads((directory / "requests.json").read_text())
    for request in requests:
        assert sha256(directory / request["response_file"]) == request["sha256"]
        payload = json.loads((directory / request["response_file"]).read_text())
        batch = [wrapper["entry"] for wrapper in payload["omim"]["entryList"]]
        assert {str(entry["mimNumber"]) for entry in batch} == set(request["requested_ids"])
        entries.extend(dict(entry, response_file=request["response_file"]) for entry in batch)
    identifiers = [str(entry["mimNumber"]) for entry in entries]
    assert len(identifiers) == len(set(identifiers))
    assert set(identifiers) == set(run["requested_ids"])
    maps, coverage = flatten_entries(entries)
    return entries, maps, coverage, requests


def load_maps(root=ROOT):
    pointer_path = root / "local-data/omim-genes/latest.json"
    if not pointer_path.exists():
        return None
    pointer = json.loads(pointer_path.read_text())
    directory = (root / pointer["run"]).resolve()
    assert directory.is_relative_to((root / "local-data/omim-genes").resolve())
    assert sha256(directory / "run.json") == pointer["run_sha256"]
    run = json.loads((directory / "run.json").read_text())
    assert sha256(root / "data/omim/omim_ids.csv") == run["seed_sha256"]
    _, targets = seed_targets(root)
    assert targets == run["requested_ids"]
    return directory, read_maps(directory)


def unique_text(values):
    return "; ".join(sorted(set(str(v) for v in values if str(v))))


def link_sources(mentions, names, maps, coverage):
    resolution = names[["omim_id", "resolved_omim_id"]]
    sources = mentions.groupby(["omim_id", "group"], sort=True).agg(
        traits=("trait", unique_text), scopes=("scope", unique_text),
        papers=("paper", unique_text), source_files=("source_file", unique_text),
        cited_genes=("cited_gene", unique_text), source_gene_sets=("source_gene_set", unique_text)
    ).reset_index().merge(resolution, on="omim_id", validate="many_to_one")
    genes = sources.merge(maps, on="resolved_omim_id", how="inner", validate="many_to_many")
    genes = genes.loc[genes["gene_symbol"].ne("")].reset_index(drop=True)
    all_ids = sources.merge(coverage, on="resolved_omim_id", how="left", validate="many_to_one")
    assert all_ids["mapping_status"].notna().all()
    return genes, all_ids


def summarize_genes(genes):
    return genes.groupby(["group", "gene_symbol"], sort=True).agg(
        cited_omim_ids=("omim_id", unique_text), resolved_omim_ids=("resolved_omim_id", unique_text),
        gene_mim_ids=("gene_mim_id", unique_text), conditions=("condition_name", unique_text),
        traits=("traits", unique_text), scopes=("scopes", unique_text),
        papers=("papers", unique_text), mapping_keys=("mapping_key", unique_text),
        response_files=("response_file", unique_text)
    ).reset_index()


def write_tables(root=ROOT):
    loaded = load_maps(root)
    if loaded is None:
        raise ValueError("Run omim_gene_maps.py --fetch first")
    directory, (_, maps, coverage, _) = loaded
    mentions = pd.read_csv(root / "data/omim/omim_id_mentions.csv", dtype=str, keep_default_na=False)
    names, _ = seed_targets(root)
    genes, all_ids = link_sources(mentions, names, maps, coverage)
    summary = summarize_genes(genes)
    for file, frame in (("gene-mappings.csv", maps), ("id-coverage.csv", all_ids),
                        ("gene-evidence.csv", genes), ("genes.csv", summary)):
        frame.to_csv(directory / file, index=False)
    for group in ("pigmentation", "hair"):
        summary.loc[summary["group"].eq(group)].to_csv(directory / f"{group}-genes.csv", index=False)
    print(json.dumps({"queried_ids": len(coverage), "unique_approved_genes": genes["gene_symbol"].nunique(),
        "genes_by_group": summary.groupby("group").size().to_dict(),
        "ids_without_approved_genes": int(coverage["approved_genes"].eq(0).sum())}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    compile_ids(ROOT)
    if args.fetch:
        fetch_maps(ROOT)
    write_tables(ROOT)


if __name__ == "__main__":
    main()
