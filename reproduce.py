#!/usr/bin/env python3
"""Rebuild OMIM and GWAS offline, pass both new gene lists to R, and export the figure.

Python 3.11+, the packages in requirements.txt, R 4.2+, and eulerr 8.3.1 are
required. All paths are relative to this file; the working directory is irrelevant.
"""

from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
GROUPS = ("pigmentation", "hair")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    require(bool(rows), f"No rows to write: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path, algorithm="sha256"):
    return hashlib.new(algorithm, path.read_bytes()).hexdigest()


def values(row, key):
    result = json.loads(row[key])
    require(isinstance(result, list) and all(isinstance(value, str) for value in result),
            f"Expected a JSON string array in {key}")
    return set(result)


def joined(items):
    return "; ".join(sorted(set(items)))


def export_gwas(output, target, counts):
    """Translate computed JSON-array columns to R's semicolon-list schema.

    Retained names come only from 08-retained-genes.csv. Evidence and mapped
    associations supply provenance; selected terms supply query labels/scopes.
    The conversion verifies membership and paper support without re-filtering.
    """
    retained = read_csv(output / "08-retained-genes.csv")
    terms = {row["efo_id"]: row for row in read_csv(output / "03-selected-terms.csv")}
    mapped_rows = read_csv(output / "06-mapped-associations.csv")
    mapped = {(row["group"], row["association_id"]): row for row in mapped_rows}
    require(len(mapped) == len(mapped_rows), "Duplicate mapped association keys")
    evidence = defaultdict(list)
    for row in read_csv(output / "07-gene-evidence.csv"):
        evidence[row["group"], row["gene_symbol"]].append(row)
    minimum = counts["filters"]["minimum_publications"]
    selection = " + ".join(counts["selection"])
    exports = {group: [] for group in GROUPS}
    seen = set()
    for gene in retained:
        group, symbol = gene["group"], gene["gene_symbol"]
        key = (group, symbol)
        require(group in GROUPS and key not in seen, "Unexpected or duplicate retained gene")
        seen.add(key)
        rows = evidence[key]
        require(bool(rows), f"No source evidence for {key}")
        association_ids = {row["association_id"] for row in rows}
        study_ids = {row["study_id"] for row in rows}
        pubmed_ids = {row["pubmed_id"] for row in rows if row["pubmed_id"]}
        for field, observed, count in (("association_ids", association_ids, "associations"),
                                       ("study_ids", study_ids, "studies"),
                                       ("pubmed_ids", pubmed_ids, "papers")):
            require(observed == values(gene, field) and len(observed) == int(gene[count]),
                    f"Retained gene differs from source evidence: {key}, {field}")
        require(gene["passes_paper_filter"].lower() == "true" and len(pubmed_ids) >= minimum,
                f"Retained gene does not satisfy the recorded paper filter: {key}")
        associations = [mapped[group, identifier] for identifier in sorted(association_ids)]
        require(all(symbol in values(row, "mapped_genes") for row in associations),
                f"Mapped association does not contain the retained name: {key}")
        for row in rows:
            association = mapped[group, row["association_id"]]
            require(all(row[field] == association[field] for field in ("study_id", "pubmed_id", "association_url"))
                    and all(values(row, field) == values(association, field)
                            for field in ("query_ids", "analysis_scopes", "source_urls")),
                    f"Gene evidence differs from its mapped association: {key}")
        combined = {field: set().union(*(values(row, field) for row in associations))
                    for field in ("query_ids", "analysis_scopes", "response_files", "source_urls")}
        query_terms = [terms[identifier] for identifier in sorted(combined["query_ids"])]
        require(all(row["group"] == group for row in query_terms), "Query belongs to another trait group")
        exports[group].append({
            "query_group": group, "gene_symbol": symbol,
            "query_ids": joined(combined["query_ids"]),
            "query_labels": joined(row["label"] for row in query_terms),
            "query_scopes": joined(row["scope"] for row in query_terms),
            "analysis_scopes": joined(combined["analysis_scopes"]),
            "associations": len(association_ids), "studies": len(study_ids),
            "publications": len(pubmed_ids),
            "reported_traits": joined(row["reported_trait"] for row in associations),
            "study_ids": joined(study_ids), "pubmed_ids": joined(pubmed_ids),
            "response_files": joined(combined["response_files"]), "selection": selection,
            "passes_publication_filter": "True", "association_ids": joined(association_ids),
            "association_urls": joined(row["association_url"] for row in associations),
            "source_urls": joined(combined["source_urls"]),
        })
    for group, rows in exports.items():
        require(len(rows) == counts["groups"][group]["retained_genes"], "Exported GWAS gene count differs")
        write_csv(target / f"gwas-{group}-genes.csv", sorted(rows, key=lambda row: row["gene_symbol"]))


def export_inputs():
    target = ROOT / "eulerr/inputs"
    target.mkdir(parents=True, exist_ok=True)
    omim = ROOT / "omim/outputs"
    gwas = ROOT / "gwas/outputs"
    counts = read_json(gwas / "figure-counts.json")
    require(set(counts["selection"]) == {"strict", "broad"}, "Unexpected GWAS selection")
    require(counts["filters"]["minimum_discovery_sample_size"] is None,
            "The R figure currently expects an unset discovery sample-size threshold")
    export_gwas(gwas, target, counts)
    manifest = []
    for group in GROUPS:
        shutil.copy2(omim / f"{group}-genes.csv", target / f"omim-{group}-genes.csv")
        for source in ("omim", "gwas"):
            path = target / f"{source}-{group}-genes.csv"
            rows = read_csv(path)
            manifest.append({"set": f"{source}_{group}", "source": source.upper(), "group": group,
                "file": path.name, "rows": len(rows), "md5": digest(path, "md5"),
                "upstream_table": f"../omim/outputs/{group}-genes.csv" if source == "omim" else
                                  "../gwas/outputs/08-retained-genes.csv",
                "selection": "cited-ID approved symbols" if source == "omim" else " + ".join(counts["selection"])})
    write_csv(target / "input-manifest.csv", manifest)
    write_csv(target / "settings.csv", [
        {"setting": "selection", "value": " + ".join(counts["selection"])},
        {"setting": "minimum_publications", "value": counts["filters"]["minimum_publications"]},
        {"setting": "minimum_discovery_sample_size", "value": "unset"},
        {"setting": "matching", "value": "Exact gene names; no alias harmonization"},
    ])
    runs = []
    sources = [("omim", "omim_" + name, path) for name, path in
               read_json(omim / "verification.json")["source_runs"].items()]
    sources += [("gwas", "gwas_" + name, path) for name, path in
                read_json(gwas / "provenance.json")["source_runs"].items()]
    for package, stage, directory in sources:
        path = ROOT / package / "inputs" / directory / "run.json"
        runs.append({"stage": stage, "run_directory": f"../{package}/inputs/{directory}",
                     "run_sha256": digest(path)})
    write_csv(target / "source-runs.csv", runs)
    return manifest


def verify_r_outputs(manifest):
    """Compare R's actual gene sets and all exclusive regions with the new inputs."""
    input_path = ROOT / "eulerr/inputs"
    output = ROOT / "eulerr/outputs"
    membership = read_csv(output / "gene-membership.csv")
    genes = {}
    for item in manifest:
        key = item["set"]
        genes[key] = {row["gene_symbol"] for row in read_csv(input_path / item["file"])}
        observed = {row["gene_symbol"] for row in membership if row[key] == "TRUE"}
        require(observed == genes[key], f"R did not use the newly generated {key} genes")
    regions = read_csv(output / "region-genes.csv")
    region_counts = read_csv(output / "overlap-regions.csv")
    supplement = read_csv(output / "euler-genes.csv")
    expected_keys = {(group, gene) for group in GROUPS for gene in
                     genes[f"omim_{group}"] | genes[f"gwas_{group}"]}
    require(len(supplement) == len(expected_keys) and
            {(row["trait_group"], row["gene_symbol"]) for row in supplement} == expected_keys,
            "Supplementary gene rows differ from the plotted genes")
    for group in GROUPS:
        omim, gwas = genes[f"omim_{group}"], genes[f"gwas_{group}"]
        expected = {"omim_only": omim - gwas, "shared": omim & gwas, "gwas_only": gwas - omim}
        for region, names in expected.items():
            actual = {row["gene_symbol"] for row in regions if row["group"] == group and row["region"] == region}
            counts = [row for row in region_counts if row["group"] == group and row["region"] == region]
            require(actual == names and len(counts) == 1 and int(counts[0]["gene_count"]) == len(names),
                    f"R region differs from newly computed genes: {group}/{region}")
            label = {"omim_only": "OMIM only", "shared": "Shared", "gwas_only": "GWAS only"}[region]
            require({row["gene_symbol"] for row in supplement if row["trait_group"] == group and
                     row["euler_region"] == label} == names,
                    f"Supplementary region differs from the plot: {group}/{region}")
    return {key: len(names) for key, names in genes.items()}


def write_file_tree():
    """List the complete local bundle, including inputs omitted by Git."""
    target = ROOT / "FILE_TREE.txt"
    target.touch()
    lines = [ROOT.name + "/"]
    excluded = {".git", ".library", ".venv", "__pycache__", ".DS_Store", "live", ".Rhistory", "Rplots.pdf"}

    def visit(directory, prefix=""):
        children = sorted((p for p in directory.iterdir() if p.name not in excluded and
                           not p.name.startswith(".env")), key=lambda p: (not p.is_dir(), p.name))
        for index, path in enumerate(children):
            last = index == len(children) - 1
            lines.append(prefix + ("└── " if last else "├── ") + path.name + ("/" if path.is_dir() else ""))
            if path.is_dir():
                visit(path, prefix + ("    " if last else "│   "))

    visit(ROOT)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    rscript = shutil.which("Rscript")
    require(rscript is not None, "Install R, then run Rscript eulerr/setup.R once")
    print("Rebuilding OMIM and GWAS from packaged snapshots (offline).", flush=True)
    for package in ("omim", "gwas"):
        subprocess.run([sys.executable, str(ROOT / package / "reproduce.py")], check=True)
    manifest = export_inputs()
    subprocess.run([rscript, "--vanilla", str(ROOT / "eulerr/reproduce.R")], check=True)
    gene_counts = verify_r_outputs(manifest)
    figures = ROOT / "figures"
    figures.mkdir(exist_ok=True)
    copied = {}
    for extension in ("pdf", "svg", "png"):
        source = ROOT / "eulerr/outputs" / f"omim-gwas-euler.{extension}"
        target = figures / source.name
        require(source.is_file() and source.stat().st_size > 0, f"Missing R figure: {source.name}")
        shutil.copy2(source, target)
        require(digest(source) == digest(target), f"Figure copy differs: {source.name}")
        copied[target.relative_to(ROOT).as_posix()] = digest(target)
    supplementary = ROOT / "supplementary-tables"
    supplementary.mkdir(exist_ok=True)
    supplementary_hashes = {}
    for source_name, target_name in (("euler-genes.csv", "euler-genes.csv"),
                                     ("euler-genes-columns.csv", "euler-genes-columns.csv"),
                                     ("euler-genes-README.md", "README.md")):
        source = ROOT / "eulerr/outputs" / source_name
        target = supplementary / target_name
        shutil.copy2(source, target)
        require(digest(source) == digest(target), f"Supplementary copy differs: {source_name}")
        supplementary_hashes[target.relative_to(ROOT).as_posix()] = digest(target)
    report = {"success": True, "mode": "offline; both source analyses recomputed",
        "analysis_api_calls": 0, "eulerr_version_required": "8.3.1", "gene_counts": gene_counts,
        "r_memberships_match_new_source_genes": True, "all_region_gene_sets_match": True,
        "inputs": [{"file": "eulerr/inputs/" + row["file"], "sha256": digest(ROOT / "eulerr/inputs" / row["file"])}
                   for row in manifest],
        "upstream_verification": {package: digest(ROOT / package / "outputs/verification.json")
                                  for package in ("omim", "gwas")},
        "figure_sha256": copied, "supplementary_sha256": supplementary_hashes,
        "driver_sha256": digest(Path(__file__))}
    (ROOT / "eulerr/outputs/pipeline-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    write_file_tree()
    print("\nVerified both newly computed gene sources, every R region, and all figure copies.")
    print(f"Final manuscript figure: {figures / 'omim-gwas-euler.pdf'}")
    print(f"Supplementary gene table: {supplementary / 'euler-genes.csv'}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        raise SystemExit(f"Reproduction failed: {error}") from None
