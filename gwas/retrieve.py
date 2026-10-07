"""Retrieve all GWAS inputs used by reproduce.py, using Python's standard library.

The five input configuration files define every query and decision. No gene
names, study IDs, selected ontology IDs, or expected counts are embedded here.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import threading
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlsplit
from urllib.request import Request, urlopen

from reproduce import (Inputs, check, sha256, save_csv, read_csv, read_pages,
                       classify_terms, read_associations, read_studies)

CONFIG_FILES = ("search-terms.json", "classification-rules.csv", "reported-trait-query.json",
                "phenotype-scope-rules.csv", "gene-filters.json")
STUDY_API = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def run_id():
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


class CatalogClient:
    """Throttle request starts across workers; retry transient HTTP errors."""
    def __init__(self):
        self.lock = threading.Lock()
        self.last_request = 0.0

    def get(self, url):
        check(urlsplit(url).scheme == "https" and urlsplit(url).netloc == "www.ebi.ac.uk",
              "Unexpected Catalog API origin")
        for attempt in range(3):
            with self.lock:
                time.sleep(max(0, .65 - (time.monotonic() - self.last_request)))
                self.last_request = time.monotonic()
            try:
                request = Request(url, headers={"Accept": "application/json", "User-Agent": "TiG-GWAS-Reproduction/1.0"})
                with urlopen(request, timeout=90) as response:
                    check(urlsplit(response.url).netloc == "www.ebi.ac.uk", "Unexpected response origin")
                    return response.read()
            except (HTTPError, URLError, TimeoutError) as error:
                if isinstance(error, HTTPError) and error.code not in (429, 500, 502, 503, 504):
                    raise
                if attempt == 2:
                    raise
                time.sleep(2 ** (attempt + 1))


def pages(client, api, parameters, collection):
    """Retrieve all pages, including an empty response; reject changing totals."""
    page_number, total, seen = 0, None, 0
    while True:
        url = api + "?" + urlencode(dict(parameters, page=page_number))
        raw = client.get(url)
        payload = json.loads(raw)
        records = payload.get("_embedded", {}).get(collection, [])
        page = payload["page"]
        check(page["number"] == page_number, "Wrong page returned")
        total = page["totalElements"] if total is None else total
        check(page["totalElements"] == total, "Catalog changed during pagination; retry this retrieval")
        seen += len(records)
        yield url, raw, records, page
        if page_number + 1 >= page["totalPages"]:
            check(seen == total, "Incomplete Catalog response")
            break
        page_number += 1


def complete_snapshot(root, directory, metadata):
    metadata["requests_sha256"] = sha256(directory / "requests.csv")
    metadata["completed_utc"] = utc_now()
    (directory / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")


def publish(root, folder, directory):
    (root / folder / "latest.json").write_text(json.dumps({
        "run": directory.relative_to(root).as_posix(),
        "run_sha256": sha256(directory / "run.json")}, indent=2) + "\n")


def fetch_terms(root, client):
    base = root / "data/gwas"
    config_path = base / "search-terms.json"
    config = json.loads(config_path.read_text())
    check(config["api"] == "https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits", "Unexpected term API")
    check(config["parameter"] == "efo_trait" and 1 <= config["page_size"] <= 500, "Invalid search configuration")
    directory = base / "runs" / run_id()
    directory.mkdir(parents=True)
    shutil.copy2(config_path, directory / "search-terms.json")
    requests, terms = [], {}
    for query in config["queries"]:
        parameters = {config["parameter"]: query["text"], "size": config["page_size"],
                      "sort": "efo_id", "direction": "asc"}
        for url, raw, records, page in pages(client, config["api"], parameters, "efo_traits"):
            file = f"response-{len(requests)+1:03d}.json"
            (directory / file).write_bytes(raw)
            requests.append(dict(query_group=query["group"], search_text=query["text"],
                page=page["number"], returned_rows=len(records), total_matches=page["totalElements"],
                url=url, response_file=file, sha256=sha256(directory/file), retrieved_utc=utc_now()))
            for term in records:
                identifier = term["efo_id"]
                value = dict(label=term["efo_trait"], uri=term["uri"])
                check(identifier not in terms or terms[identifier] == value, "Conflicting term labels")
                terms[identifier] = value
        save_csv(directory / "requests.csv", requests)
        print(f"Term search {query['text']}: {page['totalElements']} matches", flush=True)
    complete_snapshot(root, directory, dict(api=config["api"], run_id=directory.name,
        config_sha256=sha256(directory/"search-terms.json"), request_count=len(requests)))
    read_pages(Inputs(root), directory, read_csv(directory/"requests.csv"), "01_search", "efo_traits",
               ("query_group", "search_text"), "returned_rows", "total_matches")
    publish(root, "data/gwas", directory)
    return terms


def fetch_associations(root, terms, client):
    base = root / "data/gwas"
    classified = classify_terms(terms, read_csv(base/"classification-rules.csv"))
    config_path = base / "reported-trait-query.json"
    config = json.loads(config_path.read_text())
    check(config["api"] == "https://www.ebi.ac.uk/gwas/rest/api/v2/associations", "Unexpected association API")
    check(config["show_child_trait"] is False and 1 <= config["page_size"] <= 500, "Invalid association configuration")
    selected = [row for row in classified if row["group"] in config["groups"] and row["scope"] in config["scopes"]]
    check(selected, "No mapped terms selected; review the classification rules")
    selected.sort(key=lambda row: (config["groups"].index(row["group"]), config["scopes"].index(row["scope"]), row["efo_id"]))
    directory = base / "reported-labels/runs" / run_id()
    directory.mkdir(parents=True)
    shutil.copy2(config_path, directory/"query.json")
    save_csv(directory/"query-terms.csv", [{key: row[key] for key in ("efo_id", "label", "group", "scope")} for row in selected])
    requests = []
    for query in selected:
        parameters = dict(efo_id=query["efo_id"], show_child_trait="false", size=config["page_size"])
        for url, raw, records, page in pages(client, config["api"], parameters, "associations"):
            file = f"{query['efo_id']}-page-{page['number']:03d}.json"
            (directory/file).write_bytes(raw)
            requests.append(dict(query_id=query["efo_id"], query_group=query["group"], query_scope=query["scope"],
                query_label=query["label"], url=url, page=page["number"], rows=len(records),
                total_associations=page["totalElements"], response_file=file,
                sha256=sha256(directory/file), retrieved_utc=utc_now()))
        save_csv(directory/"requests.csv", requests)
        print(f"Associations {query['label']}: {page['totalElements']} records", flush=True)
    complete_snapshot(root, directory, dict(api=config["api"], query_sha256=sha256(directory/"query.json"),
        query_terms_sha256=sha256(directory/"query-terms.csv")))
    batches = read_pages(Inputs(root), directory, read_csv(directory/"requests.csv"), "04_associations", "associations",
                        ("query_id",), "rows", "total_associations")
    filters = json.loads((base/"gene-filters.json").read_text())
    _, associations = read_associations(batches, selected, read_csv(root/filters["scope_rules"]))
    publish(root, "data/gwas/reported-labels", directory)
    return directory, associations


def fetch_studies(root, association_directory, associations, client):
    seeds = sorted({row["study_id"] for row in associations})
    check(seeds, "No associations available for study retrieval")
    directory = root / "data/gwas/studies/runs" / run_id()
    directory.mkdir(parents=True)
    def one_study(identifier):
        url = STUDY_API + "/" + identifier
        study_raw = client.get(url)
        study = json.loads(study_raw)
        check(study["accession_id"] == identifier, "Wrong study accession")
        ancestry_url = study["_links"]["ancestries"]["href"].split("{")[0]
        check(ancestry_url == url + "/ancestries", "Unexpected ancestry URL")
        ancestry_raw = client.get(ancestry_url)
        ancestry = json.loads(ancestry_raw)
        check("next" not in ancestry.get("_links", {}) and "page" not in ancestry, "Paged ancestry response")
        rows = []
        for kind, body, endpoint in (("study", study_raw, url), ("ancestries", ancestry_raw, ancestry_url)):
            file = f"{identifier}-{kind}.json"
            (directory/file).write_bytes(body)
            rows.append(dict(study_id=identifier, kind=kind, url=endpoint, response_file=file,
                             sha256=sha256(directory/file), retrieved_utc=utc_now()))
        return rows
    requests = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(one_study, identifier) for identifier in seeds]
        for number, future in enumerate(as_completed(futures), 1):
            requests.extend(future.result())
            requests.sort(key=lambda row: (row["study_id"], row["kind"]))
            save_csv(directory/"requests.csv", requests)
            if number % 10 == 0 or number == len(seeds):
                print(f"Studies and sample sizes: {number}/{len(seeds)}", flush=True)
    metadata = dict(api=STUDY_API, study_ids=seeds,
        association_run=association_directory.relative_to(root).as_posix(),
        association_run_sha256=sha256(association_directory/"run.json"))
    complete_snapshot(root, directory, metadata)
    read_studies(Inputs(root), directory, metadata, read_csv(directory/"requests.csv"), associations)
    publish(root, "data/gwas/studies", directory)


def fetch_all(configuration_root, destination, client=None):
    """Build a fresh snapshot from configuration files only; preserve the pinned one."""
    configuration_root, destination = Path(configuration_root).resolve(), Path(destination).resolve()
    check(not destination.exists(), "New retrieval directory already exists; choose a different --run-directory")
    base = destination / "data/gwas"
    base.mkdir(parents=True)
    for name in CONFIG_FILES:
        shutil.copy2(configuration_root / "data/gwas" / name, base/name)
    client = CatalogClient() if client is None else client
    terms = fetch_terms(destination, client)
    directory, associations = fetch_associations(destination, terms, client)
    fetch_studies(destination, directory, associations, client)
    return destination
