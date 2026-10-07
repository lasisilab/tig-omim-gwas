# OMIM

`python reproduce.py` rebuilds gene lists from literature tables and saved responses.
`python retrieve.py` collects new responses using `OMIM_API_KEY` from the environment.

Inputs: D’Arcy Table S1, Needle Table 2, and the saved OMIM entry responses.
The code extracts cited IDs, retrieves names, resolves moved entries and reads
approved symbols from `include=geneMap` at https://api.omim.org/api/entry.

`outputs/pigmentation-genes.csv` and `outputs/hair-genes.csv` retain source IDs
and literature provenance. Other output tables record intermediate steps.

Fresh runs use `live/`; the retrieval command prints how to analyze them.
The saved OMIM payloads, outputs and reference PDF are local and ignored by Git.
