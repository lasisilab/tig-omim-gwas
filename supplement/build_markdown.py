"""Write the editable supplement from saved analysis CSVs and PubMed metadata."""
from pathlib import Path
import csv,json
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'supplement'
def read(p):
 with (ROOT/p).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def cell(v):return str(v).replace('|',r'\|').replace('\n',' ')
def table(headers,rows):
 return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(cell(x) for x in row)+' |' for row in rows])
def ids(s):return [x.strip() for x in s.split(';') if x.strip()]
refs=json.loads((OUT/'pubmed-references.json').read_text())['records']
genes=read('supplementary-tables/euler-genes.csv')
terms=read('gwas/outputs/03-selected-terms.csv')
order=[]
for row in genes:
 for pmid in ids(row['gwas_pubmed_ids']):
  if pmid not in order:order.append(pmid)
assert set(order)==set(refs)
number={p:i+3 for i,p in enumerate(order)}
lines=['# Supplemental Information','',
'## Classical foundations and genomic advances in human pigmentation and hair morphology','',
'Yemko Pryor¹, Lily Heald¹, Tina Lasisi¹˒²*','',
'¹ Department of Anthropology, University of Michigan, Ann Arbor, MI 48109, USA.  ',
'² Department of Ecology and Evolutionary Biology, University of Michigan, Ann Arbor, MI 48109, USA.','',
'*Correspondence: tlasisi@umich.edu (T. Lasisi).','',
'## Inventory of Supplemental Information','',
'- **Supplementary methods, related to Figure 3:** OMIM and GWAS Catalog retrieval, gene selection and overlap analysis.',
'- **Table S1, related to Figure 3:** Selected GWAS Catalog ontology terms and their strict/broad categories.',
'- **Table S2, related to Figure 3:** All plotted genes, their OMIM identifiers and supporting references.','',
'## Supplementary methods','',
'We extract Online Mendelian Inheritance in Man (OMIM) IDs from D’Arcy et al. [S1], Supplementary Table S1, and Needle et al. [S2], Table 2. We assign D’Arcy records and Needle hair color and graying records to pigmentation, and the remaining Needle traits to hair. We deduplicate IDs, retrieve entry names through the OMIM application programming interface (API; https://api.omim.org/api/entry), resolve moved entries, and retrieve gene mappings with `include=geneMap`. We retain the returned approved gene symbols and deduplicate them within each trait group.','',
'Next, we retrieve associations from the Genome-Wide Association Study (GWAS) Catalog v2 API (https://www.ebi.ac.uk/gwas/rest/api/v2) for the ontology IDs in Table S1. We designate core traits as strict and related traits as broad, and combine both categories without child-trait expansion. The keyword searches used to identify candidate terms are recorded in the repository. We retrieve study and ancestry records and retain mapped gene names supported by at least two distinct nonempty PubMed IDs within each trait group. We apply no sample-size threshold or additional association P-value filter. Separate papers can include overlapping cohorts; mapped genes are association annotations, not causal assignments.','',
'We match exact gene names between OMIM and the retained GWAS sets without alias harmonization. We calculate shared and source-specific sets and plot them with eulerr 8.3.1 on one area scale. Table S2 lists every plotted gene and its source identifiers and references. The saved October 3, 2026 API responses, request URLs, retrieval times and checksums preserve the analyzed snapshot.','',
'### Code and data','',
'The repository at https://github.com/lasisilab/tig-omim-gwas contains the retrieval scripts, query settings, analysis code and CSV tables. Running `python reproduce.py` regenerates the gene lists and plots from the saved responses. `supplement/build_markdown.py` generates this document and the citation-numbered Table S2 CSV. File paths below are relative to the repository root.','',
'## Table S1, related to Figure 3. Selected GWAS Catalog ontology terms','',
f'All {len(terms)} selected mapped terms are shown. Strict denotes core traits and broad denotes related traits; both categories enter the analysis. Ontology IDs are the identifiers returned by the GWAS Catalog.','',
'Source: `supplement/table-s1-selected-terms.csv`, extracted from `gwas/outputs/03-selected-terms.csv`.','',
table(['Trait group','Catalog label','Ontology ID','Scope'],[[r[k] for k in ['group','label','efo_id','scope']] for r in terms]),'',
'## Table S2, related to Figure 3. Genes represented in the Euler diagrams','',
f'All {len(genes)} gene–trait records are shown. Shared means membership in both selected sets. OMIM IDs are the literature-cited identifiers. References identify the source literature for OMIM IDs and the supporting GWAS papers. A dash indicates no identifier in the selected OMIM set. Full provenance, including PubMed and Catalog study IDs, is retained in `supplementary-tables/euler-genes.csv`.','',
'Source: `supplement/table-s2-genes-with-references.csv`, generated from `supplementary-tables/euler-genes.csv`. Citation numbers correspond to the supplemental references below.','']
rows=[]
for row in genes:
 citations=[]
 if 'Arcy' in row['omim_source_papers']:citations.append(1)
 if 'Needle' in row['omim_source_papers']:citations.append(2)
 if row['omim_cited_ids']:assert citations
 citations+= [number[p] for p in ids(row['gwas_pubmed_ids'])]
 assert citations
 rows.append([row['trait_group'],row['gene_symbol'],row['euler_region'],row['omim_cited_ids'] or '—','['+', '.join('S'+str(n) for n in sorted(set(citations)))+']'])
with (OUT/'table-s2-genes-with-references.csv').open('w',newline='') as f:
 writer=csv.writer(f);writer.writerow(['Trait group','Gene','Membership','OMIM IDs','References']);writer.writerows(rows)
with (OUT/'table-s1-selected-terms.csv').open('w',newline='') as f:
 writer=csv.writer(f);writer.writerow(['group','label','efo_id','scope']);writer.writerows([[r[k] for k in ['group','label','efo_id','scope']] for r in terms])
lines +=[table(['Trait group','Gene','Membership','OMIM IDs','References'],rows),'','## Supplemental references','',
'**S1.** D’Arcy, C. et al. (2023) Disease-gene networks of skin pigmentation disorders and reconstruction of protein-protein interaction networks. Bioengineering 10, 13. https://doi.org/10.3390/bioengineering10010013','',
'**S2.** Needle, C.D. et al. (2025) A comprehensive review of GWASs of human hair traits. J. Invest. Dermatol. 145, 2964–2972. https://doi.org/10.1016/j.jid.2025.07.004','']
for pmid in order:
 d=refs[pmid]
 author=d['author']+(' et al.' if d['multiple_authors'] else '')
 volume_pages=d['volume']+(', '+d['pages'] if d['pages'] else '')
 lines +=[f"**S{number[pmid]}.** {author} ({d['year']}) {d['title'].rstrip('.')} . {d['journal']} {volume_pages}. https://pubmed.ncbi.nlm.nih.gov/{pmid}/".replace(' .','.').replace('  ',' '),'']
path=OUT/'supplementary-information.md';path.write_text('\n'.join(lines))
assert len(rows)==473 and len(number)==49
print(f'{path}: {len(rows)} gene rows, {len(number)+2} references')
