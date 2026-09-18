# Supplementary Code

Release: CTRC VERSION 01 2026-09.

This archive contains the complete statistical analysis and plotting scripts, exact gene programs. 

## Reproduce

Use Python 3.12.10 and the pinned packages in requirements.txt, plus R 4.5.3 with survey 4.5 and survival 3.8-6. Run `python run_all.py` from the extracted directory. Rscript must be on PATH. Every script locates its inputs relative to this directory. The included results can be inspected without running any software. Runtime varies with the 2,000-resample bootstrap analyses and mixed-model fits. Ensembl coordinate annotation requests use the public GRCh37 endpoint; the retrieved coordinate snapshot is included.

The pipeline starts from the bundled patient/subtype aggregates reconstructed from all 1,254,749 deposited single cells, the NHANES raw-file merge, and the deposited bulk/spatial matrices. `prepare_raw_singlecell.py` additionally documents how to rebuild the aggregate input from GEO raw counts; it is optional because the large sparse matrix is not duplicated in this archive. Download GSE243013_NSCLC_immune_scRNA_counts.mtx.gz, GSE243013_genes.csv.gz, GSE243013_barcodes.csv.gz and GSE243013_NSCLC_immune_scRNA_metadata.csv.gz from https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/suppl/ into raw/GSE243013 before that optional step. The raw-count preparation requires enough memory for the full sparse matrix.

Data retain the terms of their original public repositories. Source clinical identifiers are deposited study labels, not names or contact details. The archive does not grant a new licence over third-party datasets.
