#!/usr/bin/env python3
"""Reproduce every count in the GWAS Catalog retrieval and retention diagram.

Run from this folder: python reproduce.py
Or run this file by its absolute path from any working directory.

Python 3.10+; standard library only. Inputs are adjacent inputs/data/gwas/:
query definitions, classification and phenotype rules, filter configuration,
and checksum-verified GWAS Catalog term, association, study and ancestry JSON.
Outputs are adjacent outputs/: diagram-counts.csv, figure-counts.json,
intermediate CSVs, gene lists, stage-input-output.csv and provenance.json.

Default: reproduce the claimed counts offline from the pinned snapshot.
python reproduce.py --fetch: retrieve new Catalog responses from the five
configuration files, then analyze them in a separate live run directory.
Both modes require only Python 3.10+; no API credentials or installed packages.
Use --inputs, --output and --run-directory to select different directories.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import parse_qs, urlsplit


def check(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def joined(values):
    return "; ".join(sorted(set(values)))


class Inputs:
    """Record every input path, checksum, stage, and original response URL."""
    def __init__(self, root):
        self.root = root
        self.files = {}

    def record(self, path, stage, expected=None, url=""):
        path = path.resolve()
        check(path.is_relative_to(self.root), "Input outside the project: " + str(path))
        digest = sha256(path)
        check(expected is None or digest == expected, "Checksum mismatch: " + str(path))
        relative = path.relative_to(self.root).as_posix()
        self.files[(stage, relative)] = {"stage": stage, "file": relative, "sha256": digest, "source_url": url}
        return path

    def json(self, path, stage, expected=None, url=""):
        return json.loads(self.record(path, stage, expected, url).read_text())

    def csv(self, path, stage, expected=None):
        return read_csv(self.record(path, stage, expected))

    def snapshot(self, folder, stage):
        pointer = self.json(self.root / folder / "latest.json", stage)
        directory = (self.root / pointer["run"]).resolve()
        check(directory.is_relative_to((self.root / folder / "runs").resolve()), "Invalid snapshot path")
        run = self.json(directory / "run.json", stage, pointer["run_sha256"])
        requests = self.csv(directory / "requests.csv", stage, run["requests_sha256"])
        return directory, run, requests


def read_pages(inputs, directory, requests, stage, collection, key_fields, rows_field, total_field):
    """Verify saved pages, completeness, and counts against the Catalog page metadata."""
    batches = defaultdict(list)
    for request in requests:
        payload = inputs.json(directory / request["response_file"], stage,
                              request["sha256"], request["url"])
        records = payload.get("_embedded", {}).get(collection, [])
        page = payload["page"]
        check(page["number"] == int(request["page"]), "Wrong page number")
        check(len(records) == int(request[rows_field]), "Wrong response row count")
        check(page["totalElements"] == int(request[total_field]), "Wrong response total")
        key = tuple(request[field] for field in key_fields)
        batches[key].append((request, page, records))
    for key, pages in batches.items():
        pages.sort(key=lambda item: item[1]["number"])
        check([p[1]["number"] for p in pages] == list(range(len(pages))), f"Missing pages: {key}")
        totals = {p[1]["totalElements"] for p in pages}
        page_totals = {p[1]["totalPages"] for p in pages}
        check(len(totals) == len(page_totals) == 1, f"Changing result totals: {key}")
        check(len(pages) == max(1, pages[0][1]["totalPages"]), f"Incomplete pagination: {key}")
        check(sum(len(p[2]) for p in pages) == pages[0][1]["totalElements"], f"Truncated results: {key}")
    return batches


def classify_terms(terms, rules):
    """Apply the first matching row of the external, ordered classification rules."""
    check(len({r["rule_id"] for r in rules}) == len(rules), "Duplicate classification rule IDs")
    result = []
    for identifier, term in sorted(terms.items()):
        matching = next((r for r in rules if re.search(r["label_regex"], term["label"], re.I)), None)
        result.append(dict(term, efo_id=identifier, **(
            {key: matching[key] for key in ("group", "scope", "rule_id")} if matching else
            {"group": "", "scope": "review", "rule_id": ""})))
    return result


def read_associations(batches, selected, scope_rules):
    """Retain query provenance, then deduplicate association IDs within each group."""
    by_id = {row["efo_id"]: row for row in selected}
    check({key[0] for key in batches} == set(by_id), "Association requests differ from selected terms")
    occurrences, unique = [], {}
    for (identifier,), pages in batches.items():
        seen = set()
        for request, _, records in pages:
            query = by_id[identifier]
            check(all(request["query_" + field] == query[field] for field in ("group", "scope", "label")),
                  "Query metadata differs from the selected term")
            for record in records:
                association_id = str(record["association_id"])
                check(association_id not in seen, "Repeated association within one ontology query")
                seen.add(association_id)
                check(identifier in {t["efo_id"] for t in record["efo_traits"]}, "Query term absent from association")
                genes = record["mapped_genes"]
                check(isinstance(genes, list) and all(isinstance(g, str) and g.strip() for g in genes),
                      "Unexpected mapped_genes value")
                row = {"group": query["group"], "query_id": identifier, "query_scope": query["scope"],
                    "association_id": association_id, "study_id": record["accession_id"],
                    "pubmed_id": str(record.get("pubmed_id") or ""), "mapped_genes": sorted(set(genes)),
                    "reported_trait": joined(record["reported_trait"]),
                    "mapped_labels": joined(t["efo_trait"] for t in record["efo_traits"]),
                    "association_url": record["_links"]["self"]["href"], "response_file": request["response_file"],
                    "source_url": request["url"], "analysis_scope": query["scope"], "scope_rule": ""}
                matching = [r for r in scope_rules if r["query_group"] == row["group"] and
                            re.search(r["pattern"], row[r["field"]], re.I)]
                check(len(matching) <= 1, "Overlapping phenotype rules")
                if matching:
                    row.update(analysis_scope=matching[0]["scope"], scope_rule=matching[0]["rule_id"])
                occurrences.append(row)
                if row["analysis_scope"] not in {"strict", "broad"}:
                    continue
                key = (row["group"], association_id)
                identity = {field: row[field] for field in ("group", "association_id", "study_id", "pubmed_id",
                    "mapped_genes", "reported_trait", "mapped_labels", "association_url")}
                if key not in unique:
                    unique[key] = dict(identity, query_ids=[], analysis_scopes=[], response_files=[], source_urls=[])
                check(all(unique[key][field] == value for field, value in identity.items()),
                      "Conflicting repeated association records")
                for field, value in (("query_ids", identifier), ("analysis_scopes", row["analysis_scope"]),
                                     ("response_files", row["response_file"]), ("source_urls", row["source_url"])):
                    unique[key][field] = sorted(set(unique[key][field]) | {value})
    return occurrences, [unique[key] for key in sorted(unique)]


def read_studies(inputs, directory, run, requests, associations):
    """Read PubMed IDs and sum initial ancestry counts per study, preserving missing N."""
    seeds = {r["study_id"] for r in associations}
    check(seeds == set(run["study_ids"]), "Study snapshot differs from selected associations")
    lookup = {(r["study_id"], r["kind"]): r for r in requests}
    check(len(lookup) == len(requests) == 2 * len(seeds), "Missing or duplicated study responses")
    studies, ancestry_rows = {}, []
    for identifier in sorted(seeds):
        payloads = {}
        for kind in ("study", "ancestries"):
            request = lookup[(identifier, kind)]
            payloads[kind] = inputs.json(directory / request["response_file"], "05_samples",
                                        request["sha256"], request["url"])
        study, ancestry = payloads["study"], payloads["ancestries"]
        check(study["accession_id"] == identifier, "Wrong study accession")
        check("next" not in ancestry.get("_links", {}) and "page" not in ancestry, "Paged ancestry response")
        records = ancestry.get("_embedded", {}).get("ancestries", [])
        for row in records:
            n = row.get("number_of_individuals")
            check(n is None or (type(n) is int and n >= 0), "Invalid ancestry sample count")
            ancestry_rows.append({"study_id": identifier, "stage": row["type"], "number_of_individuals": n,
                                 "source_url": row["_links"]["self"]["href"]})
        def count(stage):
            values = [row.get("number_of_individuals") for row in records if row["type"] == stage]
            return None if not values or None in values else sum(values)
        studies[identifier] = {"study_id": identifier, "pubmed_id": str(study.get("pubmed_id") or ""),
            "discovery_n": count("initial"), "replication_n": count("replication"),
            "source_url": lookup[(identifier, "ancestries")]["url"]}
    for row in associations:
        check(row["pubmed_id"] == studies[row["study_id"]]["pubmed_id"], "Association/study PubMed IDs disagree")
    return studies, ancestry_rows


def gene_tables(mapped_associations):
    """One evidence row per group/association/name; count distinct nonempty PubMed IDs."""
    evidence, support = [], defaultdict(lambda: defaultdict(set))
    for row in mapped_associations:
        for gene in row["mapped_genes"]:
            evidence.append({field: row[field] for field in ("group", "association_id", "study_id", "pubmed_id",
                "query_ids", "analysis_scopes", "association_url", "source_urls")} | {"gene_symbol": gene})
            sets = support[(row["group"], gene)]
            sets["association_ids"].add(row["association_id"])
            sets["study_ids"].add(row["study_id"])
            if row["pubmed_id"]:
                sets["pubmed_ids"].add(row["pubmed_id"])
    genes = [{"group": group, "gene_symbol": gene, "associations": len(sets["association_ids"]),
        "studies": len(sets["study_ids"]), "papers": len(sets["pubmed_ids"]),
        **{field: sorted(sets[field]) for field in ("association_ids", "study_ids", "pubmed_ids")}}
        for (group, gene), sets in sorted(support.items())]
    return evidence, genes


def save_csv(path, rows):
    # Array columns use JSON, preserving values instead of ambiguous comma splitting.
    fields = list(dict.fromkeys(field for row in rows for field in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v
                         for k, v in row.items()} for row in rows)


def reproduce(root, output):
    """Run every stage; write outputs only after raw-response validation succeeds."""
    inputs, tables, stages = Inputs(root), {}, []
    base = root / "data/gwas"
    def stage(name, source, operation, files):
        stages.append({"stage": name, "input": source, "operation": operation,
                       "output": "; ".join(f"{file} ({len(tables[file]):,} rows)" for file in files)})

    # 1. Read the configured keyword queries and the actual saved term-search pages.
    term_dir, term_run, term_requests = inputs.snapshot("data/gwas", "01_search")
    search = inputs.json(base / "search-terms.json", "01_search", term_run["config_sha256"])
    inputs.record(term_dir / "search-terms.json", "01_search", term_run["config_sha256"])
    batches = read_pages(inputs, term_dir, term_requests, "01_search", "efo_traits",
                         ("query_group", "search_text"), "returned_rows", "total_matches")
    queries = {(q["group"], q["text"]) for q in search["queries"]}
    check(len(queries) == len(search["queries"]) and queries == set(batches), "Search queries differ from saved requests")
    matches, terms = [], {}
    for (group, text), pages in batches.items():
        seen = set()
        for request, _, records in pages:
            check(parse_qs(urlsplit(request["url"]).query)[search["parameter"]] == [text], "Wrong search URL")
            for term in records:
                identifier = term["efo_id"]
                check(identifier not in seen, "Repeated term within one keyword query")
                seen.add(identifier)
                matches.append({"query_group": group, "search_text": text, "efo_id": identifier,
                    "label": term["efo_trait"], "uri": term["uri"], "source_url": request["url"]})
                value = {"label": term["efo_trait"], "uri": term["uri"]}
                check(identifier not in terms or terms[identifier] == value, "Conflicting ontology labels")
                terms[identifier] = value
    tables["01-search-queries.csv"] = search["queries"]
    tables["01-term-matches.csv"] = matches
    stage("1. Keyword searches", "search-terms.json + saved /efo-traits JSON pages",
          "Read every completed query and all returned pages; preserve repeated query–term matches.",
          ["01-search-queries.csv", "01-term-matches.csv"])

    # 2. Deduplicate on the returned ontology ID, not on its human-readable label.
    tables["02-mapped-terms.csv"] = [dict(efo_id=k, **v) for k, v in sorted(terms.items())]
    stage("2. Unique mapped terms", "01-term-matches.csv", "Deduplicate by efo_id.", ["02-mapped-terms.csv"])

    # 3. Apply ordered strict/broad/exclusion rules; check the actual requested seeds.
    rules = inputs.csv(base / "classification-rules.csv", "03_selection")
    classified = classify_terms(terms, rules)
    association_dir, association_run, association_requests = inputs.snapshot("data/gwas/reported-labels", "04_associations")
    query = inputs.json(base / "reported-trait-query.json", "03_selection", association_run["query_sha256"])
    inputs.record(association_dir / "query.json", "04_associations", association_run["query_sha256"])
    selected = [r for r in classified if r["group"] in query["groups"] and r["scope"] in query["scopes"]]
    recorded = inputs.csv(association_dir / "query-terms.csv", "03_selection", association_run["query_terms_sha256"])
    seed_fields = ("efo_id", "label", "group", "scope")
    check({tuple(r[k] for k in seed_fields) for r in selected} == {tuple(r[k] for k in seed_fields) for r in recorded},
          "Selected terms differ from the saved association queries; retrieve a matching snapshot")
    tables["03-classified-terms.csv"] = classified
    tables["03-selected-terms.csv"] = selected
    tables["03-unselected-terms.csv"] = [r for r in classified if r not in selected]
    stage("3. Strict + broad selection", "02-mapped-terms.csv + classification-rules.csv + reported-trait-query.json",
          "Apply the first matching classification rule; retain the configured groups and scopes.",
          ["03-classified-terms.csv", "03-selected-terms.csv", "03-unselected-terms.csv"])

    # 4. Read every selected term's association pages; deduplicate within trait group.
    config = inputs.json(base / "gene-filters.json", "04_associations")
    scope_rules = inputs.csv(root / config["scope_rules"], "04_associations")
    pages = read_pages(inputs, association_dir, association_requests, "04_associations", "associations",
                       ("query_id",), "rows", "total_associations")
    check(query["show_child_trait"] is False, "This reproduction expects queries without descendants")
    for request in association_requests:
        parameters = parse_qs(urlsplit(request["url"]).query)
        check(parameters["efo_id"] == [request["query_id"]] and parameters["show_child_trait"] == ["false"],
              "Unexpected association query parameters")
    occurrences, associations = read_associations(pages, selected, scope_rules)
    tables["04-query-association-occurrences.csv"] = occurrences
    tables["04-associations.csv"] = associations
    stage("4. Unique associations", "03-selected-terms.csv + saved /associations JSON pages + phenotype-scope-rules.csv",
          "Apply phenotype overrides; deduplicate by (group, association_id). Count distinct study IDs and nonempty PubMed IDs.",
          ["04-query-association-occurrences.csv", "04-associations.csv"])

    # 5. Read actual study/ancestry responses. The sample filter is optional.
    study_dir, study_run, study_requests = inputs.snapshot("data/gwas/studies", "05_samples")
    check(study_run["association_run_sha256"] == sha256(association_dir / "run.json"), "Study/association snapshots differ")
    studies, ancestries = read_studies(inputs, study_dir, study_run, study_requests, associations)
    inputs.record(base / "gene-filters.json", "05_samples")
    inputs.record(base / "gene-filters.json", "08_papers")
    threshold, minimum = config["minimum_discovery_sample_size"], config["minimum_publications"]
    check(threshold is None or (type(threshold) is int and threshold > 0), "Invalid discovery N threshold")
    check(type(minimum) is int and minimum > 0, "Invalid minimum paper count")
    sample_rows = [dict(row, discovery_n=studies[row["study_id"]]["discovery_n"]) for row in associations]
    eligible = [r for r in sample_rows if threshold is None or
                (r["discovery_n"] is not None and r["discovery_n"] >= threshold)]
    tables["05-study-samples.csv"] = list(studies.values())
    tables["05-ancestry-records.csv"] = ancestries
    tables["05-sample-eligible-associations.csv"] = eligible
    stage("5. Discovery sample size", "04-associations.csv + saved study/ancestry JSON + gene-filters.json",
          "Sum initial ancestry counts within each study; apply minimum_discovery_sample_size only when set.",
          ["05-study-samples.csv", "05-ancestry-records.csv", "05-sample-eligible-associations.csv"])

    # 6. Separate empty mapped_genes arrays from associations with gene names.
    mapped = [r for r in eligible if r["mapped_genes"]]
    unmapped = [r for r in eligible if not r["mapped_genes"]]
    tables["06-mapped-associations.csv"] = mapped
    tables["06-unmapped-associations.csv"] = unmapped
    stage("6. Gene mapping available", "05-sample-eligible-associations.csv",
          "Separate nonempty and empty mapped_genes arrays.",
          ["06-mapped-associations.csv", "06-unmapped-associations.csv"])

    # 7. Expand exact mapped gene strings and count names within each trait group.
    evidence, genes = gene_tables(mapped)
    tables["07-gene-evidence.csv"] = evidence
    tables["07-gene-support.csv"] = genes
    stage("7. Distinct mapped names", "06-mapped-associations.csv",
          "Expand mapped_genes; deduplicate exact names per group. No alias harmonization or extra p-value filter.",
          ["07-gene-evidence.csv", "07-gene-support.csv"])

    # 8. Count papers, not study IDs: repeated analyses from one paper count once.
    filtered = [dict(row, passes_paper_filter=row["papers"] >= minimum) for row in genes]
    retained = [r for r in filtered if r["passes_paper_filter"]]
    removed = [r for r in filtered if not r["passes_paper_filter"]]
    tables["08-retained-genes.csv"] = retained
    tables["08-removed-genes.csv"] = removed
    stage("8. Minimum-paper filter", "07-gene-evidence.csv + gene-filters.json",
          f"Retain a name with at least {minimum} distinct nonempty PubMed IDs within its trait group.",
          ["08-retained-genes.csv", "08-removed-genes.csv"])

    # Calculate every figure value from the stage objects, never from saved counts.
    counts = {"discovery": {"keyword_queries": len(queries), "unique_mapped_terms": len(terms), "selected_terms": len(selected)},
              "filters": config, "selection": query["scopes"], "groups": {}}
    for group in query["groups"]:
        subset = [r for r in associations if r["group"] == group]
        counts["groups"][group] = {"selected_terms": sum(r["group"] == group for r in selected),
            "strict_terms": sum(r["group"] == group and r["scope"] == "strict" for r in selected),
            "broad_terms": sum(r["group"] == group and r["scope"] == "broad" for r in selected),
            "strict_associations": sum("strict" in r["analysis_scopes"] for r in subset),
            "broad_associations": sum("strict" not in r["analysis_scopes"] for r in subset),
            "retrieved_associations": len(subset), "catalog_studies": len({r["study_id"] for r in subset}),
            "papers": len({r["pubmed_id"] for r in subset} - {""}),
            "sample_eligible_associations": sum(r["group"] == group for r in eligible),
            "mapped_associations": sum(r["group"] == group for r in mapped),
            "associations_without_mapped_genes": sum(r["group"] == group for r in unmapped),
            "mapped_genes": sum(r["group"] == group for r in genes),
            "retained_genes": sum(r["group"] == group for r in retained),
            "removed_genes": sum(r["group"] == group for r in removed)}
    # Reproduce the plotted branches independently using only the stage objects.
    flows = []
    def branch(section, group, source, target, value, unit):
        check(value >= 0, "Negative Sankey branch")
        flows.append(dict(section=section, group=group, source=source, target=target, value=value, unit=unit))
    for group in ("pigmentation", "hair"):
        for scope in ("strict", "broad"):
            branch("terms", group, "returned_terms", group + "_" + scope,
                   counts["groups"][group][scope + "_terms"], "ontology IDs")
    branch("terms", "", "returned_terms", "excluded_terms", len(terms)-len(selected), "ontology IDs")
    counts["discovery"]["excluded_terms"] = len(terms)-len(selected)
    for group, values in counts["groups"].items():
        for scope in ("strict", "broad"):
            value = values[scope + "_associations"]
            branch("associations", group, "retrieved", scope, value, "association IDs")
            branch("associations", group, scope, "sample_decision", value, "association IDs")
        branch("associations", group, "sample_decision", "eligible", values["sample_eligible_associations"], "association IDs")
        branch("associations", group, "sample_decision", "sample_excluded",
               values["retrieved_associations"]-values["sample_eligible_associations"], "association IDs")
        branch("associations", group, "eligible", "mapped", values["mapped_associations"], "association IDs")
        branch("associations", group, "eligible", "unmapped", values["associations_without_mapped_genes"], "association IDs")
        branch("genes", group, "distinct_names", "paper_decision", values["mapped_genes"], "gene names")
        branch("genes", group, "paper_decision", "retained", values["retained_genes"], "gene names")
        branch("genes", group, "paper_decision", "removed", values["removed_genes"], "gene names")
    tables["09-workflow-sankey-flows.csv"] = flows
    # A compact table of every number printed in the supplied diagram.
    diagram = [dict(panel="All queries", stage=label, count=counts["discovery"][key], unit=unit)
        for key, label, unit in (("keyword_queries", "Keyword queries", "queries"),
                                ("unique_mapped_terms", "Distinct mapped terms", "ontology IDs"),
                                ("selected_terms", "Selected terms", "ontology IDs"))]
    for group, values in counts["groups"].items():
        for key, label, unit in (
            ("selected_terms", "Selected mapped terms: strict + broad", "ontology IDs"),
            ("retrieved_associations", "Retrieved associations", "association IDs"),
            ("catalog_studies", "Catalog studies", "Catalog study IDs"),
            ("papers", "Papers", "distinct nonempty PubMed IDs"),
            ("mapped_associations", "With mapped gene names", "association IDs"),
            ("associations_without_mapped_genes", "No mapped genes", "association IDs"),
            ("mapped_genes", "Distinct mapped gene names", "exact gene names"),
            ("retained_genes", f"At least {minimum} papers", "exact gene names"),
            ("removed_genes", "Removed by the paper filter", "exact gene names")):
            diagram.append(dict(panel=group, stage=label, count=values[key], unit=unit))
    tables["diagram-counts.csv"] = diagram
    output.mkdir(parents=True, exist_ok=True)
    tables["stage-input-output.csv"] = stages
    tables["input-manifest.csv"] = list(inputs.files.values())
    for filename, rows in tables.items():
        save_csv(output / filename, rows)
    (output / "figure-counts.json").write_text(json.dumps(counts, indent=2) + "\n")
    provenance = {"script_sha256": sha256(Path(__file__)), "analysis_api_calls": 0,
        "source_runs": {"terms": term_dir.relative_to(root).as_posix(),
                        "associations": association_dir.relative_to(root).as_posix(),
                        "studies": study_dir.relative_to(root).as_posix()},
        "outputs": [{"file": filename, "rows": len(rows), "sha256": sha256(output / filename)}
                    for filename, rows in tables.items()]}
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return counts, stages


def print_report(counts, stages, output):
    discovery = counts["discovery"]
    print(f"{discovery['keyword_queries']} keyword queries -> {discovery['unique_mapped_terms']} distinct mapped terms"
          f" -> {discovery['selected_terms']} selected terms\n")
    groups = list(counts["groups"])
    print(f"{'Metric':<40}" + "".join(f"{group:>16}" for group in groups))
    for key in next(iter(counts["groups"].values())):
        print(f"{key:<40}" + "".join(f"{counts['groups'][group][key]:>16,}" for group in groups))
    threshold = counts["filters"]["minimum_discovery_sample_size"]
    print("\nDiscovery N threshold: " + ("unset; no exclusions for sample size" if threshold is None else str(threshold)))
    print(f"Minimum papers: {counts['filters']['minimum_publications']} distinct nonempty PubMed IDs per group/name")
    for row in stages:
        print(f"\n{row['stage']}\n  Input: {row['input']}\n  Operation: {row['operation']}\n  Output: {row['output']}")
    print(f"\nSaved tables, counts, and input checksums: {output}")
    print("Groups can overlap. Mapped names are variant annotations, not causal assignments.")


def verify_claims(counts, claims_path, output, require_exact=True):
    """Compare computed results with recorded claims; never use claims to compute."""
    expected = json.loads(claims_path.read_text())
    differences = []
    def compare(left, right, path=""):
        if isinstance(left, dict) and isinstance(right, dict):
            for key in sorted(set(left) | set(right)):
                compare(left.get(key), right.get(key), f"{path}.{key}".strip("."))
        elif left != right:
            differences.append(dict(metric=path, claimed=left, computed=right))
    compare(expected, counts)
    report = dict(matches_claimed_counts=not differences, differences=differences,
                  claims_sha256=sha256(claims_path))
    (output / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print("\nClaim verification: " + ("all values match" if not differences else f"{len(differences)} values differ"))
    if require_exact:
        check(not differences, "Computed counts differ from the pinned claims; see outputs/verification.json")
    return report


def main():
    from datetime import datetime, timezone
    bundle = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=bundle / "inputs",
                        help="Input configuration/snapshot root; defaults to adjacent inputs/")
    parser.add_argument("--output", type=Path, help="Output directory; defaults to outputs/ for pinned reproduction")
    parser.add_argument("--fetch", action="store_true", help="Retrieve terms, associations, studies and samples from the Catalog before analyzing")
    parser.add_argument("--run-directory", type=Path, help="New live-run directory; must not already exist")
    args = parser.parse_args()
    root = args.inputs.resolve()
    if args.fetch:
        from retrieve import fetch_all
        run = args.run_directory or bundle / "live" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        check(not run.resolve().exists(), "Live-run directory already exists; choose a new directory")
        root = fetch_all(root, run.resolve() / "inputs")
        output = args.output.resolve() if args.output else root.parent / "outputs"
    else:
        check(args.run_directory is None, "--run-directory requires --fetch")
        output = args.output.resolve() if args.output else bundle / "outputs"
    counts, stages = reproduce(root, output)
    print_report(counts, stages, output)
    verify_claims(counts, bundle / "inputs/claimed-counts.json", output, require_exact=not args.fetch)
    provenance_path = output / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    provenance["mode"] = "fresh API retrieval" if args.fetch else "pinned snapshot reproduction"
    provenance["retrieval_script_sha256"] = sha256(bundle / "retrieve.py")
    provenance["recorded_api_responses"] = sum(len(read_csv(root / directory / "requests.csv"))
        for directory in provenance["source_runs"].values())
    provenance["outputs"].append(dict(file="verification.json", sha256=sha256(output / "verification.json")))
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
