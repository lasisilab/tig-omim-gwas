# Genes in the Euler diagrams

`euler-genes.csv`: one row per exact gene name and trait group; all six diagram regions included.
Pigmentation: 428 names; hair: 45. The 14 shared names appear twice: 473 rows, 459 distinct names (428 + 45 - 14). Column definitions: `euler-genes-columns.csv`.

## Shared names

Final memberships below; 10 names occur in both GWAS sets and 4 in both OMIM sets. These categories overlap.

| Gene name | Pigmentation sources | Hair sources |
|---|---|---|
| ALX4 | GWAS | OMIM + GWAS |
| CIROZ | GWAS | GWAS |
| FGF5 | GWAS | GWAS |
| FRAS1 | OMIM + GWAS | OMIM + GWAS |
| HERC2 | OMIM + GWAS | GWAS |
| HOXC13 | OMIM | OMIM |
| IRF4 | OMIM + GWAS | GWAS |
| LGR4 | GWAS | OMIM |
| LINC01956 | GWAS | GWAS |
| PADI3 | OMIM + GWAS | OMIM |
| PAX3 | OMIM + GWAS | OMIM + GWAS |
| PEX14 | GWAS | OMIM |
| SLC24A4 | GWAS | GWAS |
| SLC45A2 | OMIM + GWAS | GWAS |

## Methods

- GWAS: core (strict) and related (broad) traits pooled; at least 2 distinct nonempty PubMed IDs per exact name/group. Each paper counts once. Support may pool variants/phenotypes; papers may reuse cohorts. Replication of the same variant/phenotype is not required. No sample-size cutoff.
- OMIM: cited literature tables and mapped entries; no paper-count filter. Needle et al.'s hair color/graying entries join pigmentation; other hair entries join hair (scalp, facial, and body traits).
- Overlap: exact names, without alias harmonization. Between-group overlap differs from OMIM/GWAS overlap within groups. GWAS mapped names are annotations, not causal assignments.
- Catalog scope: GCST007486 (PMID 30166351) has mixed color/morphology mappings for HERC2, IRF4, SLC24A4, and SLC45A2, contributing one of two hair-group papers for each. Group overlap does not establish causal effects on scalp-fiber geometry.

Import names/IDs as text. Multivalued fields use semicolons, also present in source labels. Blank fields mean absence from the compiled source set.

The R script verifies every region. Rebuild with the package's `reproduce.py`.
