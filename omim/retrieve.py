#!/usr/bin/env python3
"""Optional fresh OMIM API retrieval; the packaged snapshot is never overwritten.

Resolve OMIM_API_KEY in the environment, then run: python retrieve.py
Requires the same pandas/openpyxl dependencies as reproduce.py.
"""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
from omim_ids import compile_ids, fetch_names, group_table, sha256, write_json
from omim_gene_maps import fetch_maps, write_tables


def main():
    key = os.environ.get("OMIM_API_KEY", "").strip()
    if not key or key.startswith("op://"):
        raise ValueError("Supply a resolved OMIM_API_KEY through the environment")
    live = ROOT / "live" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    live.mkdir(parents=True, exist_ok=False)
    inputs = ROOT / "inputs"
    manifest = json.loads((inputs / "input-manifest.json").read_text())
    files = ["input-manifest.json", "data/omim/gene-query.json"]
    files += [row["file"] for row in manifest["files"]]
    for name in files:
        target = live / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(inputs / name, target)
    mentions, unique, _ = compile_ids(live)
    names, _ = fetch_names(unique, live)
    run = json.loads((live / "local-data/omim/latest.json").read_text())["run_directory"]
    for group in ("pigmentation", "hair"):
        group_table(mentions, names, group).to_csv(live / run / f"{group}-conditions.csv", index=False)
    fetch_maps(live)
    write_tables(live)
    write_json(live / "checksums.json", [
        {"file": p.relative_to(live).as_posix(), "sha256": sha256(p)}
        for p in sorted(live.rglob("*")) if p.is_file()
    ])
    print(f"Fresh retrieval saved in: {live}")
    output = ROOT / "outputs" / live.name
    print(f'Recalculate from it: python "{ROOT / "reproduce.py"}" --inputs "{live}" --output "{output}"')


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit(str(error)) from None
