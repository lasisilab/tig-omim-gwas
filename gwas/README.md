# GWAS Catalog

`python reproduce.py` analyzes the saved responses.
`python reproduce.py --fetch` retrieves and analyzes a new snapshot under `live/`.

`inputs/data/gwas/` contains keyword queries, classification rules, phenotype
scope overrides, filters, response JSONs and request logs.

The v2 API retrieves `/efo-traits`, `/associations`, `/studies/{accession}` and
`/studies/{accession}/ancestries` at https://www.ebi.ac.uk/gwas/rest/api/v2.
All pages are saved and checked. Strict and broad traits are combined.
Mapped names require support from two distinct PubMed IDs; no sample-size
threshold is applied. Mapped names are annotations, not causal assignments.

`outputs/08-retained-genes.csv` contains retained names and paper support.
`outputs/stage-input-output.csv` describes each step and its output tables.
