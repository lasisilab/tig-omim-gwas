# Euler plots

```sh
Rscript setup.R
Rscript reproduce.R
```

Requires eulerr 8.3.1. The R script reads the four gene CSVs in `inputs/`, checks
their manifest hashes and computes exact-name intersections. `settings.csv`
records selection and paper filtering. No APIs are called by this script.

`outputs/region-genes.csv` lists every region’s genes; `overlap-regions.csv`
gives the counts. `euler-genes.csv` joins the genes to their source evidence.
PDF, SVG and PNG plots use one area scale, with pigmentation above hair.

Run the parent `reproduce.py` to rebuild the inputs from both source analyses.
The OMIM input CSVs and outputs are local and ignored by Git.
