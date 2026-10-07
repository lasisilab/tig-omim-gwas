"""Compile cited OMIM IDs and retrieve their entry names from the OMIM API.

Run from any directory. All inputs and outputs stay inside this project.
The API key is read only from OMIM_API_KEY; a Quarto render never calls the API.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import HTTPRedirectHandler, Request, build_opener

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
API = "https://api.omim.org/api/entry"
PIGMENTATION_TRAITS = {"Hair color", "Hair graying"}
ID_PATTERN = re.compile(r"\b(\d{6})\b(?:\s*\(([^)]+)\))?")
ENTRY_TYPES = {
    "*": "Gene",
    "+": "Gene and phenotype",
    "#": "Phenotype: molecular basis known",
    "%": "Phenotype: molecular basis unknown",
    "": "Phenotype: Mendelian basis uncertain",
    "^": "Moved or removed",
}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def source_tables(project=PROJECT):
    project = Path(project)
    manifest = json.loads((project / "input-manifest.json").read_text())
    sources = [source for source in manifest["files"] if source["file"] in {
        "data/darcy2023_tableS1.xlsx", "data/needle2025_table2.csv"}]
    for source in sources:
        if sha256(project / source["file"]) != source["sha256"]:
            raise ValueError(f"Source checksum changed: {source['file']}")
    darcy_source = next(source for source in sources if source["file"].endswith(".xlsx"))
    darcy = pd.read_excel(project / darcy_source["file"],
                         sheet_name=darcy_source["worksheet"],
                         header=darcy_source["header_row"] - 1, usecols="A:G",
                         dtype=str, keep_default_na=False)
    needle = pd.read_csv(project / "data/needle2025_table2.csv", dtype=str,
                         keep_default_na=False)
    return darcy, needle, darcy_source, sources


def source_id_cells(project=PROJECT):
    """Select the two OMIM columns, keeping their original cells and provenance."""
    darcy, needle, darcy_source, sources = source_tables(project)
    darcy_cells = darcy.rename(columns={
        "Gene name": "source_gene_set", "Disease name": "source_condition",
        "Phenotype MIM number": "source_id_cell", "Phenotype class": "trait",
    })[["source_gene_set", "source_condition", "source_id_cell", "trait"]].copy()
    darcy_cells = darcy_cells.assign(
        paper=darcy_source["paper"], table=darcy_source["material"],
        source_file=darcy_source["file"], data_row=darcy.index + 1,
        source_row=darcy.index + darcy_source["header_row"] + 1,
        group="pigmentation", scope="skin pigmentation",
        source_gene_column="Gene name", source_id_column="Phenotype MIM number",
    )
    needle_source = next(source for source in sources if source["file"].endswith(".csv"))
    needle_cells = needle.rename(columns={
        "Gene Set": "source_gene_set", "Mendelian Disease Genes": "source_id_cell",
        "Trait": "trait",
    })[["source_gene_set", "source_id_cell", "trait"]].copy()
    needle_cells = needle_cells.assign(
        paper=needle_source["paper"], table=needle_source["material"],
        source_file=needle_source["file"], data_row=needle.index + 1,
        source_row=needle.index + 2, source_condition="",
        source_gene_column="Gene Set", source_id_column="Mendelian Disease Genes",
        group=needle_cells["trait"].isin(PIGMENTATION_TRAITS).map({True: "pigmentation", False: "hair"}),
        scope="broader hair traits",
    )
    needle_cells.loc[needle_cells["group"] == "pigmentation", "scope"] = "hair pigmentation"
    needle_cells.loc[needle_cells["trait"] == "Hair morphology/curliness", "scope"] = "scalp hair morphology"
    cells = pd.concat([darcy_cells, needle_cells], ignore_index=True)
    for column in ("trait", "source_gene_set", "source_condition"):
        cells[column] = cells[column].str.strip()
    return cells


def extract_ids(cells):
    """Find all six-digit IDs, then explode each source cell to one ID per row."""
    extracted = cells.assign(matches=cells["source_id_cell"].str.findall(ID_PATTERN))
    missing = extracted["matches"].str.len().eq(0)
    if (missing & ~extracted["source_id_cell"].str.strip().isin(["", "N/A"])).any():
        raise ValueError("An OMIM cell could not be parsed; inspect the source table")
    extracted = extracted.explode("matches").dropna(subset=["matches"]).copy()
    extracted["omim_id"] = extracted["matches"].str[0]
    qualifier = extracted["matches"].str[1].str.strip()
    # A gene set is not an ID-to-every-gene mapping. Use explicit qualifiers first.
    unambiguous = (extracted["source_gene_column"].eq("Gene name") |
                   ~extracted["source_gene_set"].str.contains(r"[,/\-]", regex=True))
    extracted["cited_gene"] = qualifier.where(qualifier.ne(""),
        extracted["source_gene_set"].where(unambiguous, ""))
    columns = ["omim_id", "group", "scope", "paper", "table", "source_file",
               "data_row", "source_row", "trait", "source_gene_set", "cited_gene",
               "source_condition", "source_id_cell"]
    return extracted[columns].reset_index(drop=True)


def deduplicate_ids(mentions):
    """Collapse repeated IDs while retaining all source groups and traits."""
    def unique_values(values):
        return "; ".join(sorted(set(values)))
    return mentions.groupby("omim_id", sort=True).agg(
        groups=("group", unique_values), papers=("paper", unique_values),
        traits=("trait", unique_values), scopes=("scope", unique_values),
        source_mentions=("omim_id", "size"),
    ).reset_index()


def compile_ids(project=PROJECT):
    """Write the source cells, extracted mentions, and deduplicated ID list."""
    project = Path(project)
    cells = source_id_cells(project)
    mentions = extract_ids(cells)
    unique = deduplicate_ids(mentions)
    _, _, _, sources = source_tables(project)
    output = project / "data/omim"
    output.mkdir(parents=True, exist_ok=True)
    cells.to_csv(output / "omim_source_cells.csv", index=False)
    mentions.to_csv(output / "omim_id_mentions.csv", index=False)
    unique.to_csv(output / "omim_ids.csv", index=False)
    groups = {group: set(mentions.loc[mentions["group"] == group, "omim_id"])
              for group in ("pigmentation", "hair")}
    summary = {
        "unique_ids": len(unique),
        "pigmentation_ids": len(groups["pigmentation"]), "hair_ids": len(groups["hair"]),
        "shared_ids": len(groups["pigmentation"] & groups["hair"]),
        "ids_by_source": mentions.groupby("source_file")["omim_id"].nunique().to_dict(),
        "source_mentions": len(mentions),
        "source_files": [{"file": source["file"], "sha256": source["sha256"]}
                         for source in sources],
        "classification": "D’Arcy and Needle hair color/graying: pigmentation; other Needle traits: hair",
    }
    write_json(output / "compilation.json", summary)
    return mentions, unique, summary


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError("OMIM API redirect refused")


def request_entries(ids, key, include=None):
    """Return the response and a credential-free public request URL."""
    params = {"mimNumber": ",".join(ids), "format": "json"}
    if include:
        params["include"] = include
    public_url = API + "?" + urlencode(params)
    private_url = API + "?" + urlencode(dict(params, apiKey=key))
    opener = build_opener(NoRedirect())
    for attempt in range(3):
        try:
            request = Request(private_url, headers={"Accept": "application/json",
                              "User-Agent": "gene-overlap-omim-ids/1.0"})
            with opener.open(request, timeout=30) as response:
                body = response.read()
            # The API may echo apiKey in its parameter list. Never save it.
            for secret in {key, quote(key, safe=""), json.dumps(key)[1:-1]}:
                body = body.replace(secret.encode(), b"[REDACTED]")
            payload = json.loads(body)
            if "error" in payload.get("omim", {}):
                raise ValueError("OMIM API returned an error response")
            if "entryList" not in payload.get("omim", {}):
                raise ValueError("OMIM response is missing entryList")
            return payload, public_url
        except HTTPError as error:
            code = error.code
            if code not in {429, 500, 502, 503, 504} or attempt == 2:
                raise ValueError(f"OMIM API returned HTTP {code}") from None
        except (URLError, TimeoutError, OSError):
            if attempt == 2:
                raise ValueError("OMIM API connection failed") from None
        time.sleep(2 ** attempt)
    raise ValueError("OMIM API collection incomplete")


def fetch_names(unique, project=PROJECT):
    key = os.environ.get("OMIM_API_KEY", "").strip()
    if not key or key.startswith("op://"):
        raise ValueError("Set OMIM_API_KEY to the resolved API credential before --fetch")
    project = Path(project)
    ids = unique["omim_id"].tolist()
    run_name = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run_dir = project / "local-data/omim" / run_name
    run_dir.mkdir(parents=True, exist_ok=False)
    requests, entries = [], {}
    for offset in range(0, len(ids), 20):
        batch = ids[offset:offset + 20]
        payload, url = request_entries(batch, key)
        retrieved_at = utc_now()
        response_file = f"response-{len(requests) + 1:02d}.json"
        write_json(run_dir / response_file, payload)
        returned = []
        for wrapper in payload["omim"]["entryList"]:
            entry = wrapper["entry"]
            mim = str(entry["mimNumber"])
            if mim not in batch or mim in entries:
                raise ValueError("OMIM returned an unexpected or duplicate ID")
            title = entry.get("titles", {}).get("preferredTitle", "")
            if not title:
                raise ValueError(f"OMIM returned no preferred title for {mim}")
            prefix = entry.get("prefix", "")
            if prefix not in ENTRY_TYPES:
                raise ValueError(f"OMIM returned an unknown entry prefix for {mim}")
            entries[mim] = {
                "omim_id": mim, "preferred_title": title, "prefix": prefix,
                "entry_type": ENTRY_TYPES[prefix], "status": entry.get("status", ""),
                "alternative_titles": entry.get("titles", {}).get("alternativeTitles", ""),
                "included_titles": entry.get("titles", {}).get("includedTitles", ""),
                "retrieved_at_utc": retrieved_at, "omim_url": f"https://omim.org/entry/{mim}",
            }
            returned.append(mim)
        requests.append({"url": url, "requested_ids": batch, "returned_ids": returned,
                         "retrieved_at_utc": retrieved_at, "response_file": response_file,
                         "sha256": sha256(run_dir / response_file)})
        write_json(run_dir / "requests.json", requests)
        if set(returned) != set(batch):
            raise ValueError("OMIM omitted requested IDs: " + ", ".join(sorted(set(batch) - set(returned))))
        print(f"Resolved {len(entries)}/{len(ids)} IDs", flush=True)
        if offset + 20 < len(ids):
            time.sleep(0.7)
    names = pd.DataFrame([entries[mim] for mim in ids])
    names.to_csv(run_dir / "entry-names.csv", index=False)
    metadata = {"api": API, "retrieved_at_utc": utc_now(), "unique_ids": len(names),
                "seed_sha256": sha256(project / "data/omim/omim_ids.csv"),
                "names_sha256": sha256(run_dir / "entry-names.csv"),
                "source": "Online Mendelian Inheritance in Man (OMIM), Johns Hopkins University",
                "query": "Exact cited IDs; no gene search or phenotype expansion"}
    write_json(run_dir / "run.json", metadata)
    write_json(project / "local-data/omim/latest.json",
               {"run_directory": run_dir.relative_to(project).as_posix()})
    return resolve_moved_titles(names), metadata


def resolve_moved_titles(names):
    """Keep the cited ID; supply a successor title when already retrieved."""
    names = names.copy()
    titles = dict(zip(names["omim_id"], names["preferred_title"]))
    names["condition_name"] = names["preferred_title"]
    names["resolved_omim_id"] = names["omim_id"]
    names["lookup_note"] = ""
    for index, row in names.iterrows():
        if row["prefix"] == "^":
            match = re.fullmatch(r"MOVED TO (\d{6})", row["preferred_title"])
            if match:
                target = match.group(1)
                names.loc[index, "resolved_omim_id"] = target
                names.loc[index, "lookup_note"] = f"Cited entry moved to {target}"
                if target in titles:
                    names.loc[index, "condition_name"] = titles[target]
    return names


def load_names(project=PROJECT):
    """Read the completed local lookup, refusing mismatched seed or response files."""
    project = Path(project).resolve()
    pointer = project / "local-data/omim/latest.json"
    if not pointer.exists():
        return None, None
    run_dir = (project / json.loads(pointer.read_text())["run_directory"]).resolve()
    if not run_dir.is_relative_to(project / "local-data/omim"):
        raise ValueError("OMIM run directory is outside this project")
    metadata = json.loads((run_dir / "run.json").read_text())
    if metadata["seed_sha256"] != sha256(project / "data/omim/omim_ids.csv"):
        raise ValueError("OMIM lookup uses a different ID list; rerun --fetch")
    if metadata["names_sha256"] != sha256(run_dir / "entry-names.csv"):
        raise ValueError("OMIM names checksum changed")
    for request in json.loads((run_dir / "requests.json").read_text()):
        if sha256(run_dir / request["response_file"]) != request["sha256"]:
            raise ValueError("OMIM response checksum changed")
    names = pd.read_csv(run_dir / "entry-names.csv", dtype=str, keep_default_na=False)
    return resolve_moved_titles(names), metadata


def group_table(mentions, names, group):
    """One row per ID within each source-defined trait group."""
    rows = mentions.loc[mentions["group"] == group]
    grouped = rows.groupby("omim_id", sort=True).agg(
        traits=("trait", lambda values: "; ".join(sorted(set(values)))),
        scopes=("scope", lambda values: "; ".join(sorted(set(values)))),
        papers=("paper", lambda values: "; ".join(sorted(set(values)))),
        cited_genes=("cited_gene", lambda values: "; ".join(sorted(set(values) - {""}))),
        source_gene_sets=("source_gene_set", lambda values: "; ".join(sorted(set(values)))),
    ).reset_index()
    return grouped.merge(names, on="omim_id", how="left", validate="one_to_one")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true", help="Retrieve names with OMIM_API_KEY")
    args = parser.parse_args()
    mentions, unique, summary = compile_ids()
    print(json.dumps({k: v for k, v in summary.items() if k.endswith("ids") or k == "source_mentions"}))
    names, metadata = fetch_names(unique) if args.fetch else load_names()
    if names is not None:
        run_dir = PROJECT / json.loads((PROJECT / "local-data/omim/latest.json").read_text())["run_directory"]
        for group in ("pigmentation", "hair"):
            group_table(mentions, names, group).to_csv(run_dir / f"{group}-conditions.csv", index=False)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, json.JSONDecodeError) as error:
        # Never print request objects, URLs with credentials, or raw API errors.
        raise SystemExit(str(error)) from None
