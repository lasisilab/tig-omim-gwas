# Supplemental Information

## Classical foundations and genomic advances in human pigmentation and hair morphology

Yemko Pryor¹, Lily Heald¹, Tina Lasisi¹˒²*

¹ Department of Anthropology, University of Michigan, Ann Arbor, MI 48109, USA.  
² Department of Ecology and Evolutionary Biology, University of Michigan, Ann Arbor, MI 48109, USA.

*Correspondence: tlasisi@umich.edu (T. Lasisi).

## Inventory of Supplemental Information

- **Supplementary methods, related to Figure 3:** OMIM and GWAS Catalog retrieval, gene selection and overlap analysis.
- **Table S1, related to Figure 3:** Selected GWAS Catalog ontology terms and their strict/broad categories.
- **Table S2, related to Figure 3:** All plotted genes, their OMIM identifiers and supporting references.

## Supplementary methods

We extract Online Mendelian Inheritance in Man (OMIM) IDs from D’Arcy et al. [S1], Supplementary Table S1, and Needle et al. [S2], Table 2. We assign D’Arcy records and Needle hair color and graying records to pigmentation, and the remaining Needle traits to hair. We deduplicate IDs, retrieve entry names through the OMIM application programming interface (API; https://api.omim.org/api/entry), resolve moved entries, and retrieve gene mappings with `include=geneMap`. We retain the returned approved gene symbols and deduplicate them within each trait group.

Next, we retrieve associations from the Genome-Wide Association Study (GWAS) Catalog v2 API (https://www.ebi.ac.uk/gwas/rest/api/v2) for the ontology IDs in Table S1. We designate core traits as strict and related traits as broad, and combine both categories without child-trait expansion. The keyword searches used to identify candidate terms are recorded in the repository. We retrieve study and ancestry records and retain mapped gene names supported by at least two distinct nonempty PubMed IDs within each trait group. We apply no sample-size threshold or additional association P-value filter. Separate papers can include overlapping cohorts; mapped genes are association annotations, not causal assignments.

We match exact gene names between OMIM and the retained GWAS sets without alias harmonization. We calculate shared and source-specific sets and plot them with eulerr 8.3.1 on one area scale. Table S2 lists every plotted gene and its source identifiers and references. The saved October 3, 2026 API responses, request URLs, retrieval times and checksums preserve the analyzed snapshot.

### Code and data

The repository at https://github.com/lasisilab/tig-omim-gwas contains the retrieval scripts, query settings, analysis code and CSV tables. Running `python reproduce.py` regenerates the gene lists and plots from the saved responses. `supplement/build_markdown.py` generates this document and the citation-numbered Table S2 CSV. File paths below are relative to the repository root.

## Table S1, related to Figure 3. Selected GWAS Catalog ontology terms

All 20 selected mapped terms are shown. Strict denotes core traits and broad denotes related traits; both categories enter the analysis. Ontology IDs are the identifiers returned by the GWAS Catalog.

Source: `supplement/table-s1-selected-terms.csv`, extracted from `gwas/outputs/03-selected-terms.csv`.

| Trait group | Catalog label | Ontology ID | Scope |
| --- | --- | --- | --- |
| pigmentation | hair color | EFO_0003924 | strict |
| pigmentation | eye color | EFO_0003949 | strict |
| pigmentation | freckles | EFO_0003963 | broad |
| pigmentation | suntan | EFO_0004279 | broad |
| hair | hair morphology | EFO_0005038 | strict |
| hair | synophrys measurement | EFO_0007906 | broad |
| hair | age at first facial hair | EFO_0009716 | broad |
| pigmentation | eye colour measurement | EFO_0009764 | strict |
| pigmentation | melanin measurement | EFO_0021835 | strict |
| hair | Widow's peak | HP_0000349 | broad |
| hair | Hirsutism | HP_0001007 | broad |
| pigmentation | sunburn | MONDO_0005326 | broad |
| pigmentation | skin sensitivity to sun | MONDO_0005434 | broad |
| hair | hypertrichosis | MONDO_0019280 | broad |
| pigmentation | strand of hair color | OBA_0002313 | strict |
| hair | strand of hair shape | OBA_0003493 | strict |
| pigmentation | facial pigmentation | OBA_2045282 | strict |
| hair | facial hair thickness | OBA_2050112 | broad |
| hair | coat/hair morphology trait | OBA_VT0000367 | strict |
| pigmentation | skin pigmentation | OBA_VT0002095 | strict |

## Table S2, related to Figure 3. Genes represented in the Euler diagrams

All 473 gene–trait records are shown. Shared means membership in both selected sets. OMIM IDs are the literature-cited identifiers. References identify the source literature for OMIM IDs and the supporting GWAS papers. A dash indicates no identifier in the selected OMIM set. Full provenance, including PubMed and Catalog study IDs, is retained in `supplementary-tables/euler-genes.csv`.

Source: `supplement/table-s2-genes-with-references.csv`, generated from `supplementary-tables/euler-genes.csv`. Citation numbers correspond to the supplemental references below.

| Trait group | Gene | Membership | OMIM IDs | References |
| --- | --- | --- | --- | --- |
| pigmentation | ABCB6 | OMIM only | 615402 | [S1] |
| pigmentation | ABCC9 | OMIM only | 239850 | [S1] |
| pigmentation | ACD | OMIM only | 616553 | [S1] |
| pigmentation | ACSF3 | GWAS only | — | [S3, S4, S5] |
| pigmentation | ADAM10 | OMIM only | 615537 | [S1] |
| pigmentation | ADAMTS12 | GWAS only | — | [S3, S4, S6] |
| pigmentation | ADAR | OMIM only | 127400; 615010 | [S1] |
| pigmentation | ADD3 | GWAS only | — | [S3, S7] |
| pigmentation | ADGRV1 | GWAS only | — | [S3, S7, S8] |
| pigmentation | AFG3L1P | GWAS only | — | [S4, S8, S9, S10] |
| pigmentation | AHNAK | GWAS only | — | [S3, S7, S8] |
| pigmentation | AKAP1 | GWAS only | — | [S3, S8, S11] |
| pigmentation | AKAP12 | GWAS only | — | [S3, S7, S8] |
| pigmentation | ALX4 | GWAS only | — | [S3, S8] |
| pigmentation | AMACR | GWAS only | — | [S6, S8] |
| pigmentation | ANKRD11 | GWAS only | — | [S3, S4, S5, S6, S10, S12] |
| pigmentation | AP3B1 | OMIM only | 608233 | [S1] |
| pigmentation | AP3D1 | OMIM only | 617050 | [S1] |
| pigmentation | APBA2 | GWAS only | — | [S3, S4, S6, S13] |
| pigmentation | APC2 | OMIM only | 617169 | [S1] |
| pigmentation | AREG | GWAS only | — | [S3, S8] |
| pigmentation | ARHGAP24 | GWAS only | — | [S3, S8] |
| pigmentation | ARID1B | GWAS only | — | [S3, S14] |
| pigmentation | ASIP | GWAS only | — | [S3, S15, S16] |
| pigmentation | ATM | OMIM only | 208900 | [S1] |
| pigmentation | ATP10A | GWAS only | — | [S3, S7] |
| pigmentation | ATP11A | GWAS only | — | [S3, S10, S12] |
| pigmentation | ATP1B3 | GWAS only | — | [S3, S8] |
| pigmentation | ATP2B1-AS1 | GWAS only | — | [S3, S4] |
| pigmentation | ATP5MGP4 | GWAS only | — | [S3, S8] |
| pigmentation | ATP7A | OMIM only | 309400 | [S1] |
| pigmentation | ATP7B | OMIM only | 277900 | [S1] |
| pigmentation | ATXN7L1 | GWAS only | — | [S3, S8] |
| pigmentation | AURKBP1 | GWAS only | — | [S3, S8] |
| pigmentation | AXIN2 | GWAS only | — | [S3, S8] |
| pigmentation | BANP | GWAS only | — | [S4, S6] |
| pigmentation | BCAS1 | GWAS only | — | [S3, S8, S17] |
| pigmentation | BLM | OMIM only | 210900 | [S1] |
| pigmentation | BLOC1S3 | OMIM only | 614077 | [S1] |
| pigmentation | BLOC1S5 | OMIM only | 619172 | [S1] |
| pigmentation | BLOC1S6 | OMIM only | 614171 | [S1] |
| pigmentation | BMP7 | GWAS only | — | [S3, S8, S18] |
| pigmentation | BNC1 | GWAS only | — | [S3, S8] |
| pigmentation | BNC2 | GWAS only | — | [S3, S6, S8, S10, S11, S12, S17, S19, S20, S21, S22] |
| pigmentation | BPIFA2 | GWAS only | — | [S3, S10] |
| pigmentation | BRAF | OMIM only | 613707 | [S1] |
| pigmentation | BRCA2 | OMIM only | 605724 | [S1] |
| pigmentation | BRD10 | GWAS only | — | [S3, S4, S8] |
| pigmentation | BRIP1 | OMIM only | 609054 | [S1] |
| pigmentation | BTC | GWAS only | — | [S3, S8] |
| pigmentation | C5orf67 | GWAS only | — | [S3, S7, S8] |
| pigmentation | CA5A | GWAS only | — | [S3, S4] |
| pigmentation | CBFA2T3 | GWAS only | — | [S3, S4, S5, S8] |
| pigmentation | CBS | OMIM only | 236200 | [S1] |
| pigmentation | CCND2-AS1 | GWAS only | — | [S3, S8] |
| pigmentation | CDC42BPA | GWAS only | — | [S3, S8] |
| pigmentation | CDH3 | OMIM only | 601553 | [S1] |
| pigmentation | CDK10 | GWAS only | — | [S4, S6, S8, S23] |
| pigmentation | CDK14 | GWAS only | — | [S3, S8] |
| pigmentation | CDKAL1 | GWAS only | — | [S3, S8] |
| pigmentation | CDKN2A | OMIM only | 155601; 155755; 606719 | [S1] |
| pigmentation | CDKN2B-AS1 | GWAS only | — | [S3, S7, S21] |
| pigmentation | CFAP299 | GWAS only | — | [S3, S8] |
| pigmentation | CHL1 | GWAS only | — | [S3, S8] |
| pigmentation | CHMP1A | GWAS only | — | [S3, S4, S5, S6, S8, S10] |
| pigmentation | CIROZ | GWAS only | — | [S3, S8] |
| pigmentation | CLIC5 | GWAS only | — | [S3, S8] |
| pigmentation | CNTN4 | GWAS only | — | [S6, S24] |
| pigmentation | CPNE7 | GWAS only | — | [S4, S6, S10, S25] |
| pigmentation | CPSF2 | GWAS only | — | [S3, S4, S8] |
| pigmentation | CREB5 | GWAS only | — | [S3, S8] |
| pigmentation | CTNS | OMIM only | 219800 | [S1] |
| pigmentation | CYCSP30 | GWAS only | — | [S3, S14] |
| pigmentation | CYP11A1 | OMIM only | 613743 | [S1] |
| pigmentation | CYP1B1 | GWAS only | — | [S3, S10] |
| pigmentation | DBNDD1 | GWAS only | — | [S4, S5, S8, S10, S25] |
| pigmentation | DCT | GWAS only | — | [S3, S8, S10, S17] |
| pigmentation | DDB2 | OMIM only | 278740 | [S1] |
| pigmentation | DDHD1-DT | GWAS only | — | [S3, S8] |
| pigmentation | DDIT4-AS1 | GWAS only | — | [S3, S8] |
| pigmentation | DDX3X | OMIM only | 300958 | [S1] |
| pigmentation | DEF8 | GWAS only | — | [S3, S4, S5, S6, S8, S9, S10, S12] |
| pigmentation | DISC1FP1 | GWAS only | — | [S3, S10, S14] |
| pigmentation | DKC1 | OMIM only | 305000 | [S1] |
| pigmentation | DLGAP4 | GWAS only | — | [S3, S10] |
| pigmentation | DNAI3 | GWAS only | — | [S3, S8] |
| pigmentation | DNAJB4 | GWAS only | — | [S3, S7] |
| pigmentation | DOCK8-AS1 | GWAS only | — | [S3, S7, S8] |
| pigmentation | DPEP1 | GWAS only | — | [S3, S4, S6, S25] |
| pigmentation | DSTYK | Shared | 270750 | [S1, S3, S7, S8, S17, S26] |
| pigmentation | DTNBP1 | OMIM only | 614076 | [S1] |
| pigmentation | DUSP22 | GWAS only | — | [S3, S6, S8, S10] |
| pigmentation | DUT | GWAS only | — | [S14, S27] |
| pigmentation | EDN3 | Shared | 613265 | [S1, S3, S8] |
| pigmentation | EDNRB | OMIM only | 277580; 600155; 600501 | [S1, S2] |
| pigmentation | EDNRB-AS1 | GWAS only | — | [S3, S8, S16, S17] |
| pigmentation | EGR3-AS1 | GWAS only | — | [S3, S8] |
| pigmentation | EIF6 | GWAS only | — | [S3, S19] |
| pigmentation | ENPP1 | OMIM only | 615522 | [S1] |
| pigmentation | ENTREP2 | GWAS only | — | [S3, S13] |
| pigmentation | EPG5 | OMIM only | 242840 | [S1] |
| pigmentation | EPHA4 | GWAS only | — | [S3, S8, S17] |
| pigmentation | ERCC2 | OMIM only | 278730 | [S1] |
| pigmentation | ERCC3 | OMIM only | 610651 | [S1] |
| pigmentation | ERCC4 | OMIM only | 278760 | [S1] |
| pigmentation | ERCC5 | OMIM only | 278780 | [S1] |
| pigmentation | ERCC6 | OMIM only | 600630 | [S1] |
| pigmentation | ERRFI1-DT | GWAS only | — | [S3, S8, S17] |
| pigmentation | ESCO2 | OMIM only | 268300 | [S1] |
| pigmentation | ETV1 | GWAS only | — | [S3, S8] |
| pigmentation | EXOC2 | GWAS only | — | [S3, S4, S5, S6, S8, S9, S10, S23, S25, S28] |
| pigmentation | EXOC7P1 | GWAS only | — | [S3, S8] |
| pigmentation | EYS | GWAS only | — | [S18, S29] |
| pigmentation | EZR | GWAS only | — | [S3, S8, S17] |
| pigmentation | FAM13A | GWAS only | — | [S3, S30] |
| pigmentation | FANCA | Shared | 227650 | [S1, S3, S4, S5, S6, S8, S10, S19, S31, S32] |
| pigmentation | FANCC | OMIM only | 227645 | [S1] |
| pigmentation | FANCD2 | OMIM only | 227646 | [S1] |
| pigmentation | FANCE | OMIM only | 600901 | [S1] |
| pigmentation | FANCI | OMIM only | 609053 | [S1] |
| pigmentation | FBN1 | GWAS only | — | [S14, S27] |
| pigmentation | FGF5 | GWAS only | — | [S3, S8] |
| pigmentation | FGFR3 | OMIM only | 612247 | [S1] |
| pigmentation | FIG4 | OMIM only | 611228 | [S1] |
| pigmentation | FLNA | OMIM only | 300244 | [S1] |
| pigmentation | FMN2 | GWAS only | — | [S3, S8] |
| pigmentation | FMR1 | OMIM only | 300624 | [S1] |
| pigmentation | FOSL2-AS1 | GWAS only | — | [S3, S8] |
| pigmentation | FRAS1 | Shared | 219000 | [S1, S2, S3, S8] |
| pigmentation | FREM2 | Shared | 617666 | [S1, S3, S8] |
| pigmentation | FTO | GWAS only | — | [S3, S7, S8] |
| pigmentation | FZD1 | GWAS only | — | [S3, S8] |
| pigmentation | FZD7 | GWAS only | — | [S3, S8] |
| pigmentation | GAB2 | GWAS only | — | [S3, S8, S17] |
| pigmentation | GABRG3 | GWAS only | — | [S3, S4, S8, S13] |
| pigmentation | GALNS | GWAS only | — | [S3, S4] |
| pigmentation | GAPDHP70 | GWAS only | — | [S3, S10] |
| pigmentation | GIPC2 | GWAS only | — | [S3, S6, S7] |
| pigmentation | GMDS | GWAS only | — | [S3, S14] |
| pigmentation | GNAS | OMIM only | 174800 | [S1] |
| pigmentation | GOLGA6L7 | GWAS only | — | [S6, S13] |
| pigmentation | GPNMB | OMIM only | 617920 | [S1] |
| pigmentation | GPR143 | OMIM only | 300500 | [S1] |
| pigmentation | GRIP1 | OMIM only | 617667 | [S1] |
| pigmentation | GRM5 | GWAS only | — | [S3, S6, S8, S10, S13, S25, S26, S30] |
| pigmentation | GRM5-AS1 | GWAS only | — | [S3, S10] |
| pigmentation | HDAC4 | GWAS only | — | [S3, S8] |
| pigmentation | HDAC9 | GWAS only | — | [S3, S14] |
| pigmentation | HERC2 | Shared | 227220; 615516 | [S2, S3, S4, S5, S6, S7, S8, S9, S10, S13, S14, S16, S17, S19, S23, S26, S28, S29, S30, S31, S32, S33, S34, S35, S36, S37] |
| pigmentation | HGD | OMIM only | 203500 | [S1] |
| pigmentation | HOXC13 | OMIM only | 614931 | [S2] |
| pigmentation | HPS1 | OMIM only | 203300 | [S1] |
| pigmentation | HPS3 | OMIM only | 614072 | [S1] |
| pigmentation | HPS4 | OMIM only | 614073 | [S1] |
| pigmentation | HPS5 | OMIM only | 614074 | [S1] |
| pigmentation | HPS6 | OMIM only | 614075 | [S1] |
| pigmentation | HRAS | OMIM only | 137550; 218040 | [S1] |
| pigmentation | HSPA8P5 | GWAS only | — | [S3, S8] |
| pigmentation | IDS | OMIM only | 309900 | [S1] |
| pigmentation | IDUA | OMIM only | 607014 | [S1] |
| pigmentation | IFNG | OMIM only | 613254 | [S1] |
| pigmentation | IKBKG | OMIM only | 308300 | [S1] |
| pigmentation | IL16 | GWAS only | — | [S3, S8] |
| pigmentation | INHBA | GWAS only | — | [S3, S14] |
| pigmentation | IRF4 | Shared | 611724; 621097 | [S1, S2, S3, S4, S5, S6, S8, S9, S10, S12, S16, S17, S19, S20, S23, S25, S26, S28, S30, S32, S36] |
| pigmentation | JAZF1 | GWAS only | — | [S3, S8, S36] |
| pigmentation | KCNH1 | GWAS only | — | [S3, S8] |
| pigmentation | KIAA0930 | GWAS only | — | [S3, S7, S8, S10] |
| pigmentation | KIT | OMIM only | 154800; 172800 | [S1] |
| pigmentation | KITLG | Shared | 145250; 611664; 616697; 619947 | [S1, S2, S3, S8] |
| pigmentation | KRT14 | OMIM only | 125595; 161000 | [S1] |
| pigmentation | KRT31 | GWAS only | — | [S3, S8] |
| pigmentation | KRT5 | OMIM only | 131960; 179850 | [S1] |
| pigmentation | KRT86 | GWAS only | — | [S4, S8] |
| pigmentation | KRT87P | GWAS only | — | [S4, S8] |
| pigmentation | LEF1 | GWAS only | — | [S3, S4, S8] |
| pigmentation | LGR4 | GWAS only | — | [S3, S8] |
| pigmentation | LHX2 | GWAS only | — | [S3, S8, S17] |
| pigmentation | LINC-PINT | GWAS only | — | [S3, S7, S8] |
| pigmentation | LINC00158 | GWAS only | — | [S3, S8] |
| pigmentation | LINC01117 | GWAS only | — | [S3, S8] |
| pigmentation | LINC01381 | GWAS only | — | [S3, S8] |
| pigmentation | LINC01411 | GWAS only | — | [S3, S8] |
| pigmentation | LINC01449 | GWAS only | — | [S3, S14] |
| pigmentation | LINC01491 | GWAS only | — | [S7, S14, S27, S38] |
| pigmentation | LINC01505 | GWAS only | — | [S3, S7, S8] |
| pigmentation | LINC01622 | GWAS only | — | [S3, S4] |
| pigmentation | LINC01679 | GWAS only | — | [S3, S8, S36] |
| pigmentation | LINC01877 | GWAS only | — | [S21, S27] |
| pigmentation | LINC01937 | GWAS only | — | [S3, S8] |
| pigmentation | LINC01940 | GWAS only | — | [S3, S8] |
| pigmentation | LINC01953 | GWAS only | — | [S3, S7, S24] |
| pigmentation | LINC01956 | GWAS only | — | [S3, S8] |
| pigmentation | LINC02101 | GWAS only | — | [S3, S8] |
| pigmentation | LINC02105 | GWAS only | — | [S3, S8] |
| pigmentation | LINC02120 | GWAS only | — | [S6, S27] |
| pigmentation | LINC02138 | GWAS only | — | [S5, S10] |
| pigmentation | LINC02182 | GWAS only | — | [S3, S4, S14] |
| pigmentation | LINC02225 | GWAS only | — | [S3, S27] |
| pigmentation | LINC02235 | GWAS only | — | [S3, S8] |
| pigmentation | LINC02458 | GWAS only | — | [S3, S5, S8, S9, S10, S16, S17, S22, S28] |
| pigmentation | LINC02674 | GWAS only | — | [S3, S10, S11, S12, S21, S26] |
| pigmentation | LINC02728 | GWAS only | — | [S3, S17] |
| pigmentation | LINC02751 | GWAS only | — | [S3, S8, S36] |
| pigmentation | LINC02756 | GWAS only | — | [S3, S7] |
| pigmentation | LINC02917 | GWAS only | — | [S3, S8] |
| pigmentation | LINC02932 | GWAS only | — | [S3, S8] |
| pigmentation | LINC03066 | GWAS only | — | [S3, S10] |
| pigmentation | LNCBRM | GWAS only | — | [S3, S27] |
| pigmentation | LRMDA | OMIM only | 615179 | [S1, S2] |
| pigmentation | LURAP1L-AS1 | GWAS only | — | [S3, S6, S7, S8, S10, S15, S16, S17, S26, S32, S36, S39] |
| pigmentation | LYST | OMIM only | 214500 | [S1] |
| pigmentation | LYZL1 | GWAS only | — | [S3, S14] |
| pigmentation | MAP2K1 | OMIM only | 615279 | [S1] |
| pigmentation | MAP3K1 | GWAS only | — | [S3, S7, S8] |
| pigmentation | MARK3 | GWAS only | — | [S3, S7, S8] |
| pigmentation | MC1R | Shared | 203200; 266300 | [S1, S2, S3, S4, S5, S6, S8, S12, S16, S17, S21, S26, S28, S32] |
| pigmentation | MC2R | OMIM only | 202200 | [S1] |
| pigmentation | MCM4 | OMIM only | 609981 | [S1] |
| pigmentation | MED13L | GWAS only | — | [S7, S8] |
| pigmentation | MEN1 | OMIM only | 131100 | [S1] |
| pigmentation | MFSD12 | GWAS only | — | [S13, S21, S26, S27] |
| pigmentation | MFSD12-AS1 | GWAS only | — | [S21, S26, S27] |
| pigmentation | MIR155HG | GWAS only | — | [S3, S8] |
| pigmentation | MIR4268 | GWAS only | — | [S3, S8, S17] |
| pigmentation | MIR663AHG | GWAS only | — | [S3, S7] |
| pigmentation | MITF | Shared | 103500; 193510; 608890 | [S1, S3, S8] |
| pigmentation | MLANA | GWAS only | — | [S3, S4] |
| pigmentation | MLH1 | OMIM only | 276300 | [S1] |
| pigmentation | MLPH | OMIM only | 609227 | [S1] |
| pigmentation | MRAP | OMIM only | 607398 | [S1] |
| pigmentation | MRAS | GWAS only | — | [S3, S8] |
| pigmentation | MRPL45P1 | GWAS only | — | [S3, S8] |
| pigmentation | MSH2 | OMIM only | 619096 | [S1] |
| pigmentation | MSH6 | OMIM only | 619097 | [S1] |
| pigmentation | MSI2 | GWAS only | — | [S3, S8, S11] |
| pigmentation | MSX2 | GWAS only | — | [S3, S8] |
| pigmentation | MTCO1P57 | GWAS only | — | [S14, S18] |
| pigmentation | MTCO3P41 | GWAS only | — | [S14, S18] |
| pigmentation | MTERF1 | GWAS only | — | [S3, S8] |
| pigmentation | MYADML2 | GWAS only | — | [S3, S8] |
| pigmentation | MYEF2 | GWAS only | — | [S3, S40] |
| pigmentation | MYH7B | GWAS only | — | [S3, S10] |
| pigmentation | MYH8 | OMIM only | 608837 | [S1] |
| pigmentation | MYO1B-AS1 | GWAS only | — | [S3, S8] |
| pigmentation | MYO5A | OMIM only | 214450 | [S1] |
| pigmentation | Metazoa_SRP | GWAS only | — | [S3, S7, S14] |
| pigmentation | NBN | OMIM only | 251260 | [S1] |
| pigmentation | NCOA6 | GWAS only | — | [S6, S10] |
| pigmentation | NEK6 | GWAS only | — | [S3, S8, S17] |
| pigmentation | NF1 | OMIM only | 162200; 193520 | [S1] |
| pigmentation | NFIA | GWAS only | — | [S3, S8] |
| pigmentation | NFIX | OMIM only | 602535; 614753 | [S1] |
| pigmentation | NHP2 | OMIM only | 613987 | [S1] |
| pigmentation | NIPAL3 | GWAS only | — | [S3, S7] |
| pigmentation | NNT | OMIM only | 614736 | [S1] |
| pigmentation | NOP10 | OMIM only | 224230 | [S1] |
| pigmentation | NOX4 | GWAS only | — | [S3, S6, S10] |
| pigmentation | NPLOC4 | GWAS only | — | [S3, S8, S29, S32, S34, S36] |
| pigmentation | NR0B1 | OMIM only | 300200 | [S1] |
| pigmentation | NRAS | OMIM only | 137550 | [S1] |
| pigmentation | NSD1 | OMIM only | 117550 | [S1] |
| pigmentation | NYAP2 | GWAS only | — | [S3, S14] |
| pigmentation | OBI1-AS1 | GWAS only | — | [S3, S5] |
| pigmentation | OCA2 | Shared | 203200; 227220 | [S1, S2, S3, S4, S5, S6, S8, S9, S10, S12, S13, S16, S18, S22, S23, S26, S27, S32, S36, S37, S38, S41, S42] |
| pigmentation | OVOL1 | GWAS only | — | [S3, S8] |
| pigmentation | PABPN1L | GWAS only | — | [S4, S5] |
| pigmentation | PADI3 | Shared | 191480 | [S2, S3, S4, S8] |
| pigmentation | PAH | OMIM only | 261600 | [S1] |
| pigmentation | PALB2 | OMIM only | 610832 | [S1] |
| pigmentation | PARN | OMIM only | 616353 | [S1] |
| pigmentation | PAX3 | Shared | 148820; 193500 | [S1, S3, S8] |
| pigmentation | PCNT | OMIM only | 210720 | [S1] |
| pigmentation | PCNX1 | GWAS only | — | [S3, S30] |
| pigmentation | PDCD6IPP2 | GWAS only | — | [S3, S23] |
| pigmentation | PDE4B | GWAS only | — | [S3, S10] |
| pigmentation | PDE4D | GWAS only | — | [S3, S36] |
| pigmentation | PEBP4 | GWAS only | — | [S3, S8] |
| pigmentation | PEX14 | GWAS only | — | [S3, S8] |
| pigmentation | PIEZO1 | GWAS only | — | [S3, S4] |
| pigmentation | PIGPP3 | GWAS only | — | [S3, S29] |
| pigmentation | PIGU | GWAS only | — | [S7, S8, S10] |
| pigmentation | PKHD1 | GWAS only | — | [S3, S8] |
| pigmentation | PLXNB2 | GWAS only | — | [S3, S7, S8] |
| pigmentation | PMS2 | OMIM only | 619101 | [S1] |
| pigmentation | POFUT1 | OMIM only | 615327 | [S1] |
| pigmentation | POGLUT1 | OMIM only | 615696 | [S1] |
| pigmentation | POLA1 | OMIM only | 301220 | [S1] |
| pigmentation | POLH | OMIM only | 278750 | [S1] |
| pigmentation | POMC | Shared | 609734 | [S1, S3, S8] |
| pigmentation | PPARGC1B | GWAS only | — | [S3, S10, S11, S12, S21, S25, S43] |
| pigmentation | PPFIBP2 | GWAS only | — | [S3, S7, S8, S17] |
| pigmentation | PPM1A | GWAS only | — | [S3, S8, S17] |
| pigmentation | PRDM7 | GWAS only | — | [S3, S4, S5, S10] |
| pigmentation | PRELID1P5 | GWAS only | — | [S3, S8] |
| pigmentation | PRKAR1A | OMIM only | 160980 | [S1] |
| pigmentation | PSENEN | OMIM only | 613736 | [S1] |
| pigmentation | PTCH1 | GWAS only | — | [S3, S14] |
| pigmentation | PTEN | OMIM only | 158350 | [S1] |
| pigmentation | PTPN11 | OMIM only | 151100 | [S1] |
| pigmentation | PTPN14 | GWAS only | — | [S3, S7] |
| pigmentation | RAB11FIP2 | GWAS only | — | [S3, S10, S11, S12, S21, S26] |
| pigmentation | RAB27A | OMIM only | 607624 | [S1] |
| pigmentation | RAD51B | GWAS only | — | [S3, S8] |
| pigmentation | RAF1 | OMIM only | 611553; 611554 | [S1] |
| pigmentation | RALY | GWAS only | — | [S3, S8, S10, S12, S17, S20, S31] |
| pigmentation | RASA3 | GWAS only | — | [S3, S7] |
| pigmentation | RECQL4 | OMIM only | 268400 | [S1] |
| pigmentation | RIT1 | OMIM only | 615355 | [S1] |
| pigmentation | RN7SL720P | GWAS only | — | [S3, S6, S8, S10, S12, S17, S20] |
| pigmentation | RNF166 | GWAS only | — | [S3, S5, S10] |
| pigmentation | RNF5P1 | GWAS only | — | [S3, S8] |
| pigmentation | RNU1-117P | GWAS only | — | [S3, S5, S8, S9, S10, S16, S17, S22, S28] |
| pigmentation | RNU4-64P | GWAS only | — | [S14, S24] |
| pigmentation | RNU6-1313P | GWAS only | — | [S3, S8] |
| pigmentation | RNU6-211P | GWAS only | — | [S3, S8] |
| pigmentation | RNU6-366P | GWAS only | — | [S3, S4, S5, S8, S10, S16, S17, S19, S23, S28, S29, S32, S36] |
| pigmentation | RNU6-440P | GWAS only | — | [S3, S8] |
| pigmentation | RNU6-921P | GWAS only | — | [S3, S7, S8] |
| pigmentation | RNU6-929P | GWAS only | — | [S3, S8, S18] |
| pigmentation | RNU7-197P | GWAS only | — | [S3, S8] |
| pigmentation | RPL13 | GWAS only | — | [S5, S6, S32] |
| pigmentation | RPL18AP17 | GWAS only | — | [S14, S29] |
| pigmentation | RPL23AP68 | GWAS only | — | [S3, S14] |
| pigmentation | RPL3 | GWAS only | — | [S3, S8] |
| pigmentation | RPL34-DT | GWAS only | — | [S3, S8] |
| pigmentation | RPL6P14 | GWAS only | — | [S3, S8] |
| pigmentation | RPL7AP79 | GWAS only | — | [S3, S8] |
| pigmentation | RPS2P1 | GWAS only | — | [S3, S15, S16] |
| pigmentation | RSPO2 | GWAS only | — | [S3, S8] |
| pigmentation | RTEL1 | OMIM only | 615190 | [S1] |
| pigmentation | SAMD9 | OMIM only | 617053 | [S1] |
| pigmentation | SAMHD1 | GWAS only | — | [S3, S10] |
| pigmentation | SAMMSON | GWAS only | — | [S4, S7] |
| pigmentation | SCARB1 | GWAS only | — | [S3, S22] |
| pigmentation | SEMA6D | GWAS only | — | [S6, S14, S30] |
| pigmentation | SGIP1 | GWAS only | — | [S3, S10] |
| pigmentation | SGK1 | GWAS only | — | [S3, S8] |
| pigmentation | SGPL1 | OMIM only | 617575 | [S1] |
| pigmentation | SH3GL3 | GWAS only | — | [S3, S8] |
| pigmentation | SHC4 | GWAS only | — | [S3, S14] |
| pigmentation | SHOC2 | OMIM only | 607721 | [S1] |
| pigmentation | SIK1 | GWAS only | — | [S3, S8, S36] |
| pigmentation | SLA2 | GWAS only | — | [S3, S10] |
| pigmentation | SLC12A1 | GWAS only | — | [S14, S27] |
| pigmentation | SLC12A9 | GWAS only | — | [S3, S4, S8, S36] |
| pigmentation | SLC17A5 | OMIM only | 269920 | [S1] |
| pigmentation | SLC24A4 | GWAS only | — | [S3, S4, S5, S7, S8, S9, S10, S16, S17, S19, S23, S28, S29, S32, S36] |
| pigmentation | SLC24A5 | Shared | 113750 | [S1, S3, S4, S7, S13, S14, S26, S27, S30, S36, S38, S40, S41] |
| pigmentation | SLC29A3 | OMIM only | 602782 | [S1] |
| pigmentation | SLC38A2 | GWAS only | — | [S3, S8] |
| pigmentation | SLC45A2 | Shared | 227240; 606574 | [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S12, S13, S17, S19, S23, S25, S26, S29, S30, S31, S32, S36, S38] |
| pigmentation | SLC49A4 | GWAS only | — | [S3, S8] |
| pigmentation | SLC6A17 | GWAS only | — | [S3, S12] |
| pigmentation | SMARCAL1 | OMIM only | 242900 | [S1] |
| pigmentation | SMIM38 | GWAS only | — | [S5, S10] |
| pigmentation | SMO | OMIM only | 601707 | [S1] |
| pigmentation | SNX16 | GWAS only | — | [S3, S8] |
| pigmentation | SORBS2 | GWAS only | — | [S14, S24] |
| pigmentation | SOX10 | OMIM only | 609136; 611584; 613266 | [S1] |
| pigmentation | SOX5 | GWAS only | — | [S3, S36] |
| pigmentation | SOX6 | GWAS only | — | [S3, S4, S8] |
| pigmentation | SP2-AS1 | GWAS only | — | [S3, S8, S17] |
| pigmentation | SP6 | GWAS only | — | [S3, S8, S17] |
| pigmentation | SPAG16-DT | GWAS only | — | [S3, S7, S24] |
| pigmentation | SPATA2L | GWAS only | — | [S4, S6, S8] |
| pigmentation | SPATA33 | GWAS only | — | [S4, S5, S8, S20, S32] |
| pigmentation | SPG7 | GWAS only | — | [S4, S5, S6, S7, S32] |
| pigmentation | SPIRE2 | GWAS only | — | [S3, S4, S6, S8, S10, S12] |
| pigmentation | SPRED1 | OMIM only | 611431 | [S1] |
| pigmentation | SRC | GWAS only | — | [S3, S10] |
| pigmentation | ST3GAL5 | OMIM only | 609056 | [S1] |
| pigmentation | STAR | OMIM only | 201710 | [S1] |
| pigmentation | STK11 | OMIM only | 175200 | [S1] |
| pigmentation | SYNE2 | GWAS only | — | [S3, S7, S8] |
| pigmentation | TAC4-AS1 | GWAS only | — | [S3, S8] |
| pigmentation | TACC1 | GWAS only | — | [S3, S8] |
| pigmentation | TBCD | GWAS only | — | [S3, S8] |
| pigmentation | TCF25 | GWAS only | — | [S4, S8, S12] |
| pigmentation | TERC | OMIM only | 127550 | [S1] |
| pigmentation | TERT | OMIM only | 613989 | [S1] |
| pigmentation | TESHL | GWAS only | — | [S3, S6] |
| pigmentation | TGDS | GWAS only | — | [S3, S8, S10, S17] |
| pigmentation | THNSL2 | GWAS only | — | [S3, S8] |
| pigmentation | TINF2 | OMIM only | 613990 | [S1] |
| pigmentation | TIPARP-AS1 | GWAS only | — | [S3, S12] |
| pigmentation | TMCC2 | GWAS only | — | [S8, S17, S29] |
| pigmentation | TMEM163 | GWAS only | — | [S3, S7, S8] |
| pigmentation | TMEM258 | GWAS only | — | [S3, S8] |
| pigmentation | TMTC3 | GWAS only | — | [S3, S8, S14] |
| pigmentation | TP63 | Shared | 103285 | [S1, S3, S14] |
| pigmentation | TPCN2 | Shared | 612267 | [S1, S3, S4, S5, S7, S8, S9, S10, S15, S17, S19, S29] |
| pigmentation | TPM3P2 | GWAS only | — | [S3, S29] |
| pigmentation | TPRG1 | GWAS only | — | [S3, S14, S30] |
| pigmentation | TRPC4AP | GWAS only | — | [S3, S10] |
| pigmentation | TRPS1 | GWAS only | — | [S3, S10, S12] |
| pigmentation | TSC1 | OMIM only | 191100 | [S1] |
| pigmentation | TSC2 | OMIM only | 613254 | [S1] |
| pigmentation | TSPAN10 | GWAS only | — | [S3, S8, S32, S36] |
| pigmentation | TUBB3 | GWAS only | — | [S3, S4, S8, S12] |
| pigmentation | TUBB4BP4 | GWAS only | — | [S3, S7] |
| pigmentation | TWIST2 | GWAS only | — | [S3, S8] |
| pigmentation | TYR | Shared | 203100; 601800; 606952 | [S1, S2, S3, S6, S7, S8, S10, S12, S13, S16, S17, S19, S25, S26, S28, S36] |
| pigmentation | TYRP1 | Shared | 203290 | [S1, S2, S6, S7] |
| pigmentation | U6 | GWAS only | — | [S10, S14] |
| pigmentation | USP9X | OMIM only | 300968 | [S1] |
| pigmentation | UVSSA | OMIM only | 614640 | [S1] |
| pigmentation | VGLL4 | GWAS only | — | [S3, S8] |
| pigmentation | VPS9D1 | GWAS only | — | [S4, S6, S10, S12] |
| pigmentation | VPS9D1-AS1 | GWAS only | — | [S4, S6] |
| pigmentation | WHAMMP2 | GWAS only | — | [S3, S23] |
| pigmentation | WNT7B | GWAS only | — | [S8, S29] |
| pigmentation | WRAP53 | OMIM only | 613988 | [S1] |
| pigmentation | WRN | OMIM only | 277700 | [S1] |
| pigmentation | XPA | OMIM only | 278700 | [S1] |
| pigmentation | XPC | OMIM only | 278720 | [S1] |
| pigmentation | Y_RNA | GWAS only | — | [S3, S6, S8, S14] |
| pigmentation | ZBTB38 | GWAS only | — | [S3, S8, S17] |
| pigmentation | ZC3H18 | GWAS only | — | [S3, S4] |
| pigmentation | ZFP36L1 | GWAS only | — | [S3, S7, S8] |
| pigmentation | ZFPM1 | GWAS only | — | [S3, S4] |
| pigmentation | ZFYVE16 | GWAS only | — | [S3, S8] |
| pigmentation | ZMIZ1 | GWAS only | — | [S3, S8] |
| pigmentation | ZMPSTE24 | OMIM only | 608612 | [S1] |
| pigmentation | ZNF276 | GWAS only | — | [S4, S6, S12] |
| pigmentation | ZNF778 | GWAS only | — | [S3, S4] |
| pigmentation | ZNF831 | GWAS only | — | [S3, S8] |
| pigmentation | ZNG1A | GWAS only | — | [S3, S7, S8] |
| hair | ALX4 | Shared | 609597; 613451; 615529 | [S2, S11, S44, S45] |
| hair | BCL2 | GWAS only | — | [S11, S45] |
| hair | CCDST | GWAS only | — | [S46, S47] |
| hair | CIR1P1 | GWAS only | — | [S44, S48] |
| hair | CIROZ | GWAS only | — | [S11, S30, S45] |
| hair | EDAR | Shared | 129490; 224900 | [S2, S11, S30, S37, S44, S49] |
| hair | FGF5 | GWAS only | — | [S44, S45] |
| hair | FOXL2 | OMIM only | 110100; 608996 | [S2] |
| hair | FOXL2NB | GWAS only | — | [S30, S50] |
| hair | FOXO1 | OMIM only | 268220 | [S2] |
| hair | FOXP2 | OMIM only | 602081 | [S2] |
| hair | FRAS1 | Shared | 219000 | [S2, S46, S51] |
| hair | GATA3 | OMIM only | 146255 | [S2] |
| hair | GMDS-DT | GWAS only | — | [S43, S44] |
| hair | HAND1 | GWAS only | — | [S44, S48] |
| hair | HERC2 | GWAS only | — | [S32, S44] |
| hair | HOXC13 | OMIM only | 614931 | [S2] |
| hair | IRF4 | GWAS only | — | [S32, S44] |
| hair | KRT8P28 | GWAS only | — | [S46, S47] |
| hair | LGR4 | OMIM only | 615311; 619613 | [S2] |
| hair | LIMS1 | GWAS only | — | [S11, S50] |
| hair | LINC00708 | GWAS only | — | [S30, S51] |
| hair | LINC01432 | GWAS only | — | [S32, S44] |
| hair | LINC01956 | GWAS only | — | [S44, S47] |
| hair | LINC02230 | GWAS only | — | [S30, S44, S50] |
| hair | MBL2 | OMIM only | 614372 | [S2] |
| hair | PADI3 | OMIM only | 191480 | [S2] |
| hair | PAX3 | Shared | 122880; 148820; 193500; 268220 | [S2, S30, S44] |
| hair | PAX7 | OMIM only | 268220 | [S2] |
| hair | PEX14 | OMIM only | 614887 | [S2] |
| hair | PRR23A | GWAS only | — | [S30, S50] |
| hair | RPTN | GWAS only | — | [S47, S51] |
| hair | S100A11 | GWAS only | — | [S46, S47] |
| hair | SLC24A4 | GWAS only | — | [S32, S44] |
| hair | SLC39A12 | GWAS only | — | [S44, S45] |
| hair | SLC45A1 | OMIM only | 617532 | [S2] |
| hair | SLC45A2 | GWAS only | — | [S32, S44] |
| hair | SOX2 | OMIM only | 206900 | [S2] |
| hair | SOX2-OT | GWAS only | — | [S11, S44, S50] |
| hair | SPTLC1P4 | GWAS only | — | [S46, S47] |
| hair | TBX15 | Shared | 260660 | [S2, S11, S45] |
| hair | TCHH | Shared | 617252 | [S2, S19, S30, S46, S47, S51] |
| hair | THEM4 | GWAS only | — | [S46, S47] |
| hair | TMEM174 | GWAS only | — | [S30, S44, S50] |
| hair | WNT10A | Shared | 150400; 224750; 257980 | [S2, S19, S46] |

## Supplemental references

**S1.** D’Arcy, C. et al. (2023) Disease-gene networks of skin pigmentation disorders and reconstruction of protein-protein interaction networks. Bioengineering 10, 13. https://doi.org/10.3390/bioengineering10010013

**S2.** Needle, C.D. et al. (2025) A comprehensive review of GWASs of human hair traits. J. Invest. Dermatol. 145, 2964–2972. https://doi.org/10.1016/j.jid.2025.07.004

**S3.** Kichaev, G. et al. (2019) Leveraging Polygenic Functional Enrichment to Improve GWAS Power. Am J Hum Genet 104, 65-75. https://pubmed.ncbi.nlm.nih.gov/30595370/

**S4.** Jiang, L. et al. (2021) A generalized linear mixed model association tool for biobank-scale data. Nat Genet 53, 1616-1621. https://pubmed.ncbi.nlm.nih.gov/34737426/

**S5.** Lona-Durazo, F. et al. (2021) A large Canadian cohort provides insights into the genetic architecture of human hair colour. Commun Biol 4, 1253. https://pubmed.ncbi.nlm.nih.gov/34737440/

**S6.** Farré, X. et al. (2023) Skin Phototype and Disease: A Comprehensive Genetic Approach to Pigmentary Traits Pleiotropy Using PRS in the GCAT Cohort. Genes (Basel) 14, 149. https://pubmed.ncbi.nlm.nih.gov/36672889/

**S7.** Landi, M.T. et al. (2020) Genome-wide association meta-analyses combining multiple risk phenotypes provide insights into the genetic architecture of cutaneous melanoma susceptibility. Nat Genet 52, 494-504. https://pubmed.ncbi.nlm.nih.gov/32341527/

**S8.** Morgan, M.D. et al. (2018) Genome-wide study of hair colour in UK Biobank explains most of the SNP heritability. Nat Commun 9, 5271. https://pubmed.ncbi.nlm.nih.gov/30531825/

**S9.** Lin, B.D. et al. (2015) Heritability and Genome-Wide Association Studies for Hair Color in a Dutch Twin Family Based Sample. Genes (Basel) 6, 559-76. https://pubmed.ncbi.nlm.nih.gov/26184321/

**S10.** Visconti, A. et al. (2018) Genome-wide association study in 176,678 Europeans reveals genetic loci for tanning response to sun exposure. Nat Commun 9, 1684. https://pubmed.ncbi.nlm.nih.gov/29739929/

**S11.** Endo, C. et al. (2018) Genome-wide association study in Japanese females identifies fifteen novel skin-related trait associations. Sci Rep 8, 8974. https://pubmed.ncbi.nlm.nih.gov/29895819/

**S12.** Vollenbrock, C.E. et al. (2022) Genome-wide association study identifies novel loci associated with skin autofluorescence in individuals without diabetes. BMC Genomics 23, 840. https://pubmed.ncbi.nlm.nih.gov/36536295/

**S13.** Lona-Durazo, F. et al. (2019) Meta-analysis of GWA studies provides new insights on the genetic architecture of skin pigmentation in recently admixed populations. BMC Genet 20, 59. https://pubmed.ncbi.nlm.nih.gov/31315583/

**S14.** Jonnalagadda, M. et al. (2019) A Genome-Wide Association Study of Skin and Iris Pigmentation among Individuals of South Asian Ancestry. Genome Biol Evol 11, 1066-1076. https://pubmed.ncbi.nlm.nih.gov/30895295/

**S15.** Sulem, P. et al. (2008) Two newly identified genetic determinants of pigmentation in Europeans. Nat Genet 40, 835-7. https://pubmed.ncbi.nlm.nih.gov/18488028/

**S16.** Zhang, M. et al. (2013) Genome-wide association studies identify several new loci associated with pigmentation traits and skin cancer risk in European Americans. Hum Mol Genet 22, 2948-59. https://pubmed.ncbi.nlm.nih.gov/23548203/

**S17.** Hysi, P.G. et al. (2018) Genome-wide association meta-analysis of individuals of European ancestry identifies new loci explaining a substantial fraction of hair color variation and heritability. Nat Genet 50, 652-656. https://pubmed.ncbi.nlm.nih.gov/29662168/

**S18.** Rawofi, L. et al. (2017) Genome-wide association study of pigmentary traits (skin and iris color) in individuals of East Asian ancestry. PeerJ 5, e3951. https://pubmed.ncbi.nlm.nih.gov/29109912/

**S19.** Eriksson, N. et al. (2010) Web-based, participant-driven studies yield novel genetic associations for common traits. PLoS Genet 6, e1000993. https://pubmed.ncbi.nlm.nih.gov/20585627/

**S20.** Jacobs, L.C. et al. (2015) A Genome-Wide Association Study Identifies the Skin Color Genes IRF4, MC1R, ASIP, and BNC2 Influencing Facial Pigmented Spots. J Invest Dermatol 135, 1735-1742. https://pubmed.ncbi.nlm.nih.gov/25705849/

**S21.** Shin, J.G. et al. (2021) GWAS Analysis of 17,019 Korean Women Identifies the Variants Associated with Facial Pigmented Spots. J Invest Dermatol 141, 555-562. https://pubmed.ncbi.nlm.nih.gov/32835660/

**S22.** Seo, J.Y. et al. (2022) GWAS Identifies Multiple Genetic Loci for Skin Color in Korean Women. J Invest Dermatol 142, 1077-1084. https://pubmed.ncbi.nlm.nih.gov/34648798/

**S23.** Han, J. et al. (2008) A genome-wide association study identifies novel alleles associated with hair color and skin pigmentation. PLoS Genet 4, e1000074. https://pubmed.ncbi.nlm.nih.gov/18483556/

**S24.** Yoo, H.Y. et al. (2022) A Genome-Wide Association Study and Machine-Learning Algorithm Analysis on the Prediction of Facial Phenotypes by Genotypes in Korean Women. Clin Cosmet Investig Dermatol 15, 433-445. https://pubmed.ncbi.nlm.nih.gov/35313536/

**S25.** Nan, H. et al. (2009) Genome-wide association study of tanning phenotype in a population of European ancestry. J Invest Dermatol 129, 2250-7. https://pubmed.ncbi.nlm.nih.gov/19340012/

**S26.** Adhikari, K. et al. (2019) A GWAS in Latin Americans highlights the convergent evolution of lighter skin pigmentation in Eurasia. Nat Commun 10, 358. https://pubmed.ncbi.nlm.nih.gov/30664655/

**S27.** Hansen, M.E.B. et al. (2026) Anthropometric and cardio-metabolic trait variation and genetic associations in sub-Saharan Africa. Am J Hum Genet 113, 2001-2020. https://pubmed.ncbi.nlm.nih.gov/42586059/

**S28.** Sulem, P. et al. (2007) Genetic determinants of hair, eye and skin pigmentation in Europeans. Nat Genet 39, 1443-52. https://pubmed.ncbi.nlm.nih.gov/17952075/

**S29.** Xie, Z. et al. (2024) iGWAS: Image-based genome-wide association of self-supervised deep phenotyping of retina fundus images. PLoS Genet 20, e1011273. https://pubmed.ncbi.nlm.nih.gov/38728357/

**S30.** Adhikari, K. et al. (2016) A genome-wide association scan in admixed Latin Americans identifies loci influencing facial and scalp hair features. Nat Commun 7, 10815. https://pubmed.ncbi.nlm.nih.gov/26926045/

**S31.** Liu, F. et al. (2015) Genetics of skin color variation in Europeans: genome-wide association studies with functional follow-up. Hum Genet 134, 823-35. https://pubmed.ncbi.nlm.nih.gov/25963972/

**S32.** Galván-Femenía, I. et al. (2018) Multitrait genome association analysis identifies new susceptibility genes for human anthropometric variation in the GCAT cohort. J Med Genet 55, 765-778. https://pubmed.ncbi.nlm.nih.gov/30166351/

**S33.** Kayser, M. et al. (2008) Three genome-wide association studies and a linkage analysis identify HERC2 as a human iris color gene. Am J Hum Genet 82, 411-23. https://pubmed.ncbi.nlm.nih.gov/18252221/

**S34.** Liu, F. et al. (2010) Digital quantification of human eye color highlights genetic association of three new loci. PLoS Genet 6, e1000934. https://pubmed.ncbi.nlm.nih.gov/20463881/

**S35.** Candille, S.I. et al. (2012) Genome-wide association studies of quantitatively measured skin, hair, and eye pigmentation in four European populations. PLoS One 7, e48294. https://pubmed.ncbi.nlm.nih.gov/23118974/

**S36.** Simcoe, M. et al. (2021) Genome-wide association study in almost 195,000 individuals identifies 50 previously unidentified genetic loci for eye color. Sci Adv 7, eabd1239. https://pubmed.ncbi.nlm.nih.gov/33692100/

**S37.** Inoue, Y. et al. (2021) Search for genetic loci involved in the constitution and skin type of a Japanese women using a genome-wide association study. Exp Dermatol 30, 1787-1793. https://pubmed.ncbi.nlm.nih.gov/34265127/

**S38.** Batai, K. et al. (2021) Genetic loci associated with skin pigmentation in African Americans and their effects on vitamin D deficiency. PLoS Genet 17, e1009319. https://pubmed.ncbi.nlm.nih.gov/33600456/

**S39.** Kenny, E.E. et al. (2012) Melanesian blond hair is caused by an amino acid change in TYRP1. Science 336, 554. https://pubmed.ncbi.nlm.nih.gov/22556244/

**S40.** Martin, A.R. et al. (2017) An Unexpectedly Complex Architecture for Skin Pigmentation in Africans. Cell 171, 1340-1353.e14. https://pubmed.ncbi.nlm.nih.gov/29195075/

**S41.** Wang, F. et al. (2022) A Genome-Wide Scan on Individual Typology Angle Found Variants at SLC24A2 Associated with Skin Color Variation in Chinese Populations. J Invest Dermatol 142, 1223-1227.e14. https://pubmed.ncbi.nlm.nih.gov/34570997/

**S42.** Moon, H. et al. (2025) Comprehensive Profiling of Genetic and Nongenetic Factors that Influence Skin Traits in Asian Women from 4 Countries. J Invest Dermatol 145, 2272-2280.e10. https://pubmed.ncbi.nlm.nih.gov/40010489/

**S43.** Wang, P. et al. (2022) Novel genetic associations with five aesthetic facial traits: A genome-wide association study in the Chinese population. Front Genet 13, 967684. https://pubmed.ncbi.nlm.nih.gov/36035146/

**S44.** Pickrell, J.K. et al. (2016) Detection and interpretation of shared genetic influences on 42 human traits. Nat Genet 48, 709-17. https://pubmed.ncbi.nlm.nih.gov/27182965/

**S45.** Takala, J.H. et al. (2026) Hirsutism beyond PCOS: Genome-wide evidence for genetic factors. J Invest Dermatol 146, 2769-2778.e10. https://pubmed.ncbi.nlm.nih.gov/42190866/

**S46.** Medland, S.E. et al. (2009) Common variants in the trichohyalin gene are associated with straight hair in Europeans. Am J Hum Genet 85, 750-5. https://pubmed.ncbi.nlm.nih.gov/19896111/

**S47.** Ho, Y.Y.W. et al. (2020) Comparison of Genome-Wide Association Scans for Quantitative and Observational Measures of Human Hair Curvature. Twin Res Hum Genet 23, 271-277. https://pubmed.ncbi.nlm.nih.gov/33190678/

**S48.** Luo, J. et al. (2023) GWASs Identify Genetic Loci Associated with Human Scalp Hair Whorl Direction. J Invest Dermatol 143, 2065-2068.e10. https://pubmed.ncbi.nlm.nih.gov/37565938/

**S49.** Wu, S. et al. (2016) Genome-wide scans reveal variants at EDAR predominantly affecting hair straightness in Han Chinese and Uyghur populations. Hum Genet 135, 1279-1286. https://pubmed.ncbi.nlm.nih.gov/27487801/

**S50.** Wu, S. et al. (2018) Genome-wide association studies and CRISPR/Cas9-mediated gene editing identify regulatory variants influencing eyebrow thickness in humans. PLoS Genet 14, e1007640. https://pubmed.ncbi.nlm.nih.gov/30248107/

**S51.** Liu, F. et al. (2018) Meta-analysis of genome-wide association studies identifies 8 novel loci involved in shape variation of human head hair. Hum Mol Genet 27, 559-575. https://pubmed.ncbi.nlm.nih.gov/29220522/
