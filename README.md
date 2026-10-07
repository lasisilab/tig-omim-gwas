# Pigmentation and hair genes in OMIM and the GWAS Catalog

Code and data for Figure 3 and its supplementary tables.

## Reproduce

Requires Python 3.11+, R 4.2+, and eulerr 8.3.1.

```sh
python -m pip install -r requirements.txt
Rscript eulerr/setup.R
python reproduce.py
```

The offline run recalculates gene lists and overlaps from saved responses.
GWAS selection uses strict + broad traits and at least two distinct PubMed IDs.
Gene names are matched exactly; aliases are not harmonized.

- `omim/`: literature tables, OMIM retrieval and gene mapping.
- `gwas/`: query settings, Catalog retrieval and gene filtering.
- `eulerr/`: four gene lists and R plotting code.
- `figures/`: PDF, SVG and PNG plots.
- `supplementary-tables/euler-genes.csv`: all plotted genes and source evidence.
- `supplement/`: Word supplement, compact CSV tables and their generator.

## Fresh API retrieval

```sh
python gwas/reproduce.py --fetch
python omim/retrieve.py
```

OMIM requires `OMIM_API_KEY` in the environment. New responses go into `live/`
folders; they do not replace the pinned snapshot or automatically update the plot.

## Supplement

```sh
python -m pip install -r supplement/requirements.txt
python supplement/build.py
```

Tables are generated from the analysis CSVs. The Word file contains methods,
searches, selection rules and every plotted gene.

The local folder includes OMIM snapshots and a reference PDF excluded from Git.
A clone needs those inputs restored to repeat the complete offline analysis.
The R-only plot also requires the ignored OMIM gene CSVs. No credentials are stored.
