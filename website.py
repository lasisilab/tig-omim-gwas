"""Refresh website counts and evidence from the distributed analysis CSV."""
from pathlib import Path
import csv
import html
import re
import shutil

ROOT = Path(__file__).resolve().parent
for name in ["table-s1-selected-terms.csv", "supplementary-information.md"]:
    source = ROOT / "supplement" / name
    if source.exists():
        (ROOT / "downloads").mkdir(exist_ok=True)
        shutil.copy2(source, ROOT / "downloads" / name)
rows = list(csv.DictReader((ROOT / "supplementary-tables/euler-genes.csv").open()))
if not rows:
    raise ValueError("Gene evidence is empty")
keys = [(r["trait_group"], r["gene_symbol"]) for r in rows]
if len(keys) != len(set(keys)):
    raise ValueError("Duplicate gene–trait rows")
groups = {t: {r["gene_symbol"] for r in rows if r["trait_group"] == t}
          for t in ["pigmentation", "hair"]}
shared = sorted(groups["pigmentation"] & groups["hair"])

def replace_region(path, name, content):
    source = path.read_text()
    pattern = rf"<!-- {name}:start -->.*?<!-- {name}:end -->"
    updated, count = re.subn(pattern, lambda _: f"<!-- {name}:start -->\n\n{content}\n\n<!-- {name}:end -->", source, flags=re.S)
    if count != 1:
        raise ValueError(f"Missing or repeated {name} region in {path.name}")
    path.write_text(updated)

counts = ["| Trait group | OMIM only | Both sources | GWAS only | Total |", "|---|---:|---:|---:|---:|"]
for trait, names in groups.items():
    numbers = [sum(r["trait_group"] == trait and r["euler_region"] == region for r in rows)
               for region in ["OMIM only", "Shared", "GWAS only"]]
    counts.append(f"| {trait.capitalize()} | " + " | ".join(map(str, [*numbers, len(names)])) + " |")
counts.append(f"\nThe groups share {len(shared)} gene names: **{len(set().union(*groups.values()))} distinct names** across **{len(rows)} gene–trait rows**. Between-group overlap differs from OMIM/GWAS overlap within each group.")
replace_region(ROOT / "index.qmd", "counts", "\n".join(counts))
replace_region(ROOT / "genes.qmd", "shared", "## Names shared between trait groups\n\n" + ", ".join(shared) + ".")

def links(value, base):
    ids = [v.strip() for v in value.split(";") if v.strip()]
    if any(not v.isdigit() for v in ids):
        raise ValueError(f"Unexpected identifier: {value}")
    return ", ".join(f'<a href="{base}{v}">{v}</a>' for v in ids) or "—"

body = []
for row in rows:
    trait, membership = html.escape(row["trait_group"]), html.escape(row["euler_region"])
    cells = [html.escape(row["gene_symbol"]), trait, membership,
             links(row["omim_cited_ids"], "https://omim.org/entry/"),
             links(row["gwas_pubmed_ids"], "https://pubmed.ncbi.nlm.nih.gov/")]
    body.append(f'<tr data-trait="{trait}" data-membership="{membership}">' + "".join(f"<td>{cell}</td>" for cell in cells) + "</tr>")
p = ROOT / "genes.qmd"
s, count = re.subn(r"<tbody>.*?</tbody>", lambda _: "<tbody>" + "\n".join(body) + "</tbody>", p.read_text(), flags=re.S)
if count != 1:
    raise ValueError("Evidence table not found")
s = re.sub(r'(<p id="gene-count" aria-live="polite">).*?(</p>)', lambda m: m[1] + f"{len(rows)} rows" + m[2], s)
p.write_text(s)
print(f"Website evidence: {len(rows)} rows, {len(shared)} shared names.")
