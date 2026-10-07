# Genes in the Euler diagrams

`euler-genes.csv` contains 473 gene-trait rows and 459 distinct gene names.
Each row represents one gene in one trait group. Genes found in both groups occur twice.
All six regions are included. `euler-genes-columns.csv` defines every column.

GWAS uses strict and broad traits and requires at least 2 distinct PubMed IDs per gene and trait group.
OMIM has no paper filter. A blank source field means the gene is absent from that selected source set; it does not establish absence from the database.
Multiple identifiers and URLs are separated by semicolons. Text labels are preserved from the source tables; semicolons can also occur within labels.
Gene names and identifiers should be imported as text. GWAS mapped genes are association annotations, not causal assignments.

The R script joins the four validated source gene lists to the exact plotted memberships, verifies every region, and writes this table.
Run the package's root reproduce.py to rebuild the upstream OMIM/GWAS inputs and regenerate the table.
