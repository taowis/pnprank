# data/

Only this file and `truth/` are in git. Everything else sits under the data root set in
`config/paths.yaml` (about 19 GB at present: 13 GB `raw/`, 4.4 GB `interim/`, 0.9 GB
`external/`, 31 MB `processed/`) and is either downloaded or rebuilt by the Snakemake
workflow. This file lists every untracked file, where it came from and how to get it back.
The order in which the pipeline uses them is in `docs/workflow.md`; licences and
verification URLs are in `docs/data_sources.md`.

To rebuild from an empty data root:

```
pnprank build-candidates          # truth_sources/ snapshots (see the note below)
snakemake -c5 fetch_external      # reference files under external/
snakemake -c5                     # downloads the atlases as needed, then runs the pipeline
```

Every download is checked against a recorded checksum and the workflow stops on a
mismatch. A source that has since changed upstream will therefore fail loudly rather
than give silently different results.

## truth/ (in git)

| File | Contents |
|---|---|
| `neuropathy_genes.tsv` | The truth set: 242 genes with tier, discovery year, PMID and source URL |
| `candidates.tsv` | The 249 candidates from PanelApp and ClinGen that the truth set was curated from |
| `curation.json` | Per-gene evidence record behind each discovery year |
| `review.tsv` | Manual review decisions on earlier-paper screen hits |
| `sources.md` | How the truth set was built and which candidates were excluded |

## raw/: single-cell atlases

One directory per atlas, named by atlas ID. All files were downloaded on 2026-10-07.
The SHA-256 of each file is in `config/atlases.yaml`; rule `download_atlas_file`
(`workflow/scripts/download_atlas.py`) downloads a missing file and checks it against
that value. The CELLxGENE file names are dataset version IDs. If CELLxGENE retires a
version, find the dataset on the collection page and compare the SHA-256 before using it.
Only cells from healthy, control or naive donors are used (see `standardize.py`).

| Atlas | File | Size | Accession | Download URL |
|---|---|---|---|---|
| heming2025_human_sural | `GSE285983_RAW.tar` | 1.4 GB | GEO GSE285983 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285983/suppl/GSE285983_RAW.tar |
| heming2025_human_sural | `GSE285983_metadata_all.csv.gz` | 22 MB | GEO GSE285983 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE285nnn/GSE285983/suppl/GSE285983_metadata_all.csv.gz |
| yang2022_human_tg | `GSE197289_snRNA-seq_human_barcode_meta.csv.gz` | 842 kB | GEO GSE197289 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197289/suppl/GSE197289_snRNA-seq_human_barcode_meta.csv.gz |
| yang2022_human_tg | `GSE197289_snRNA-seq_human_raw_counts.RDS.gz` | 182 MB | GEO GSE197289 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197289/suppl/GSE197289_snRNA-seq_human_raw_counts.RDS.gz |
| yang2022_mouse_tg | `GSE197289_snRNA-seq_mouse_barcode_meta.csv.gz` | 2 MB | GEO GSE197289 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197289/suppl/GSE197289_snRNA-seq_mouse_barcode_meta.csv.gz |
| yang2022_mouse_tg | `GSE197289_snRNA-seq_mouse_raw_counts.RDS.gz` | 400 MB | GEO GSE197289 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197289/suppl/GSE197289_snRNA-seq_mouse_raw_counts.RDS.gz |
| jung2023_human_drg | `GSE201586_humanallcelltypes.rds.gz` | 4.2 GB | GEO GSE201586 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201586/suppl/GSE201586_humanallcelltypes.rds.gz |
| jung2023_human_drg | `GSE201586_series_matrix.txt.gz` | 4 kB | GEO GSE201586 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201586/matrix/GSE201586_series_matrix.txt.gz |
| sparc465_human_sympathetic | `primary/filtered_feature_bc_matrix/matrix.mtx.gz` | 2.2 GB | SPARC 465, DOI 10.26275/oqxx-si9h | https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com/465/files/primary/filtered_feature_bc_matrix/matrix.mtx.gz |
| sparc465_human_sympathetic | `primary/filtered_feature_bc_matrix/features.tsv.gz` | 3 MB | SPARC 465, DOI 10.26275/oqxx-si9h | https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com/465/files/primary/filtered_feature_bc_matrix/features.tsv.gz |
| sparc465_human_sympathetic | `primary/filtered_feature_bc_matrix/barcodes.tsv.gz` | 423 kB | SPARC 465, DOI 10.26275/oqxx-si9h | https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com/465/files/primary/filtered_feature_bc_matrix/barcodes.tsv.gz |
| sparc465_human_sympathetic | `primary/meta.data.csv` | 12 MB | SPARC 465, DOI 10.26275/oqxx-si9h | https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com/465/files/primary/meta.data.csv |
| yadav2023_human_spinal_cord | `GSE190442_aggregated_counts_postqc.csv.gz` | 129 MB | GEO GSE190442 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190442/suppl/GSE190442_aggregated_counts_postqc.csv.gz |
| yadav2023_human_spinal_cord | `GSE190442_aggregated_metadata_postqc.csv.gz` | 1 MB | GEO GSE190442 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190442/suppl/GSE190442_aggregated_metadata_postqc.csv.gz |
| kedlian2024_human_muscle | `5fd93c55-8047-48fb-905c-6ef2f7767483.h5ad` | 1.2 GB | CELLxGENE 2d40e6a7-f2fd-49ba-9db9-6b97e4c6dad5 | https://datasets.cellxgene.cziscience.com/5fd93c55-8047-48fb-905c-6ef2f7767483.h5ad |
| yim2022_mouse_nerves (not used, see decisions) | `GSE182098_sciatic-nerve_single-nuclei-atlas-16Dec2020.RDS.gz` | 3.1 GB | GEO GSE182098 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE182nnn/GSE182098/suppl/GSE182098_sciatic-nerve_single-nuclei-atlas-16Dec2020.RDS.gz |
| stefanovska2025_mouse_sciatic | `GSE280857_RAW.tar` | 174 MB | GEO GSE280857 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE280nnn/GSE280857/suppl/GSE280857_RAW.tar |
| stefanovska2025_mouse_sciatic | `all_final.h5ad` | 410 MB | Zenodo 14900774 | https://zenodo.org/records/14900774/files/all_final.h5ad |
| krauter2025_mouse_drg | `GSE287551_aggregated_filtered_feature_bc_matrix.h5` | 257 MB | GEO GSE287551 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287551/suppl/GSE287551_aggregated_filtered_feature_bc_matrix.h5 |
| krauter2025_mouse_drg | `GSE287551_metadata.tsv.gz` | 360 kB | GEO GSE287551 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287551/suppl/GSE287551_metadata.tsv.gz |
| alkaslasi2021_mouse_mn | `8b1fc940-a0e4-463e-a2e4-7007974fb736.h5ad` | 242 MB | CELLxGENE 101d14fc-754d-42ad-883d-d810b00b06e3 | https://datasets.cellxgene.cziscience.com/8b1fc940-a0e4-463e-a2e4-7007974fb736.h5ad |
| alkaslasi2021_mouse_mn | `GSE167597_all_AlkaslasiPiccus_annotations.csv.gz` | 334 kB | GEO GSE167597 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167597/suppl/GSE167597_all_AlkaslasiPiccus_annotations.csv.gz |
| ham2025_mouse_nmj | `b69ec11e-bab4-4e12-8c52-0028fc3c2ed3.h5ad` | 304 MB | CELLxGENE ef9f0fae-c12c-4f1a-bb67-63346c2ae0f1 | https://datasets.cellxgene.cziscience.com/b69ec11e-bab4-4e12-8c52-0028fc3c2ed3.h5ad |


## external/: reference and annotation files

| Path | Contents | Source | Version, date retrieved | How to obtain |
|---|---|---|---|---|
| `gnomad/gnomad.v4.1.constraint_metrics.tsv` | Gene constraint metrics (pLI, LOEUF, missense Z) | https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv | v4.1, 2026-10-07 | `fetch_external` (MD5 check) |
| `gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz` | Median TPM per gene and tissue | https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz | v10, 2026-10-07 | `fetch_external` (MD5 check) |
| `annotation/Mus_musculus.GRCm39.116.gtf.gz` | Mouse gene annotation | https://ftp.ensembl.org/pub/release-116/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz | Ensembl 116, 2026-10-07 | `fetch_external` (BSD sum check) |
| `annotation/mouse_genes_e116.tsv` | Mouse gene ID, name and biotype, from the GTF | derived | Ensembl 116 | `fetch_external` |
| `orthology/mus_musculus/Compara.116.protein_default.homologies.tsv.gz` | Ensembl Compara homologies, mouse dump | https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/mus_musculus/Compara.116.protein_default.homologies.tsv.gz | Ensembl 116, 2026-10-07 | `fetch_external` (MD5 check) |
| `orthology/human_mouse_one2one_e116.tsv` | 16,335 human-mouse one-to-one orthologue pairs | derived from the Compara dump | Ensembl 116 | `fetch_external` |
| `loci/cytoBand.hg38.txt.gz` | UCSC hg38 cytoband coordinates | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz | host file dated 2022-10-29, retrieved 2026-10-08 | `fetch_external` (SHA-256 check) |
| `loci/genes_<locus>.tsv` (7 files) | Protein-coding genes overlapping each linkage interval in `config/loci.tsv` | https://rest.ensembl.org/overlap/region/human/ and `/lookup/symbol/` | Ensembl 116 (stored in each file), 2026-10-08 | written by `locus_case_study` on first run, then reused |
| `truth_sources/hgnc_complete_set.txt` | HGNC gene symbols, aliases and Ensembl IDs | https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt | snapshot 2026-10-06 | `pnprank build-candidates` |
| `truth_sources/panelapp_85.json`, `panelapp_846.json`, `panelapp_165.json` | Genomics England PanelApp panels (hereditary neuropathy; neuropathy or pain GMS; HSP) | https://panelapp.genomicsengland.co.uk/api/v1/panels/{85,846,165}/ | v1.513, v8.37, v1.316; 2026-10-06 | `pnprank build-candidates` |
| `truth_sources/clingen_gene_validity.csv` | ClinGen gene-disease validity curations | https://search.clinicalgenome.org/kb/gene-validity/download | file created 2026-10-06 | `pnprank build-candidates` |
| `truth_sources/manifest.json` | URL, date and SHA-256 of the four snapshots above | | | `pnprank build-candidates` |
| `comparators/goa/*.gz` (6 files) | EBI GOA human annotation releases 25, 58, 92, 128, 163, 197 | https://ftp.ebi.ac.uk/pub/databases/GO/goa/old/HUMAN/ | release dates in `config/comparators.yaml`; 2026-10-08 | `comparators_fetch` |
| `comparators/go/<release>/*.obo` (6 files) | GO ontology releases 2004-12-01 to 2019-12-09 | https://release.geneontology.org/ | 2026-10-08 | `comparators_fetch` |
| `comparators/biogrid/*.zip` (5 files) | BioGRID ORGANISM releases 2.0.36, 3.1.72, 3.2.108, 3.4.144, 3.5.180 | https://biogrid-downloads.nyc3.digitaloceanspaces.com/BioGRID/Release-Archive/ | 2026-10-08 | `comparators_fetch` |
| `comparators/uniprot/2026_03/sec_ac.txt`, `HUMAN_9606_idmapping.dat.gz` | UniProt secondary accessions and human ID mapping (identifier maps for the old GOA files) | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/ | UniProt 2026_03, 2026-10-08 | `comparators_fetch` |
| `comparators/manifest.json` | Size, MD5 and SHA-256 of every comparator file | | | `comparators_fetch` |
| `celegans/` | Cached JSON responses from HGNC REST, Ensembl REST (pan-taxonomic Compara) and the Alliance of Genome Resources API, plus `top200_all_orthology.tsv` | URL patterns in `docs/data_sources.md` | Ensembl REST 116, Alliance 9.1.0; 2026-10-08 | `celegans_orthologues` (live queries; reuses the cache) |
| `atlas_search/` | Stage 2 search records: CELLxGENE Discover dataset listing (`cxg_datasets.json`) and the per-search notes summarised in `results/tables/atlas_candidates.tsv` | CELLxGENE Discover curation API; GEO, Zenodo, SPARC, ArrayExpress | 2026-10-07 | not rebuilt by the workflow; kept as the search record |

The truth-source snapshots cannot be fetched again as they were: HGNC, PanelApp and
ClinGen serve only their current versions, and `pnprank build-candidates` also rewrites
`truth/candidates.tsv`. To reproduce the analysis, only `hgnc_complete_set.txt` is needed
(for symbol-to-ID mapping); a newer snapshot can map a handful of renamed symbols
differently. The SHA-256 of the files used is in `truth_sources/manifest.json`. The
same applies to the UniProt identifier maps, which come from `current_release`.

Three files under `orthology/` are left over from setting up orthology and are not used:
`Compara.116.protein_default.homologies.tsv.gz` (the homo_sapiens dump, which in
release 116 has no mouse rows), its `MD5SUM`, and `ens_human_mouse_raw.tsv` (an empty
result from an Ensembl BioMart query that failed). They can be deleted.

## interim/: per-atlas intermediates

Written by `standardize_atlas` (`workflow/scripts/standardize.py`), one directory per atlas:

| File | Contents |
|---|---|
| `cells.h5ad` | Standardised cells: raw counts, species-native Ensembl IDs, author cell type, donor, keep flag |
| `pseudobulk_sums.tsv.gz` | Summed counts per gene and (cell type, donor) group; mouse genes mapped to human one-to-one orthologues |
| `pseudobulk_groups.tsv` | Cells per group, used for the minimum-cells and minimum-donors filters |
| `counts/`, `gex/`, `h5/` | Format conversions for some atlases (Matrix Market from the R scripts, gene-expression rows from the SPARC multiome matrix, per-sample h5 files from the Heming tar) |

`interim/yim2022_mouse_nerves/counts/` is from the excluded Yim atlas.
`interim/comparators/` holds the per-cutoff comparator scores written by
`comparators_backtest`.

## processed/

`features_static.tsv.gz` is the gene-by-feature matrix (19,229 protein-coding genes;
atlas specificity, GTEx and gnomAD features) written by `stage4_model`. It contains no
label-derived features; those are recomputed inside each cutoff and fold from the genes
known at that point and are not stored.

## Other untracked outputs

`results/backtest/scores/` holds the out-of-fold scores for every cutoff, written by
`stage5_backtest` and read by `comparators_backtest` and `locus_case_study`. `logs/`
holds run logs and the leakage-test record. Both are rebuilt by the workflow.
