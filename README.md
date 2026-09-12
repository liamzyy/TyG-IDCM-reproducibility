# TyG-IDCM: GSE141910 analysis materials

This repository contains the analysis code and complete moderated-t ranked gene list for the GSE141910 secondary transcriptomic analysis in the TyG-IDCM study. It does not contain the clinical dataset or clinical analysis code.

## Files

- `21_gse141910_analysis.R`: workflow for differential expression (DCM versus non-failing controls) and pre-ranked Hallmark/KEGG gene set enrichment analysis.
- `GSE141910_ranked_list.csv`: complete post-filter ranked gene universe, with 18,861 unique gene symbols and their limma moderated t-statistics, sorted in descending order. Positive values indicate higher expression in DCM. This is not a significance-filtered gene subset or a leading-edge list.

## Running the analysis

Use R with the packages `GEOquery`, `limma`, `edgeR`, `fgsea`, `msigdbr`, and `dplyr` installed. From the repository root, create `GSE141910_input/` and place these public input files inside it:

```text
GSE141910_input/
  GSE141910_raw_counts_GRCh38.p13_NCBI.tsv.gz
  Human.GRCh38.p13.annot.tsv.gz
```

The source dataset is available from [GEO: GSE141910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910). The script retrieves phenotype metadata from GEO; it does not download the expression and annotation files automatically.

Run:

```sh
Rscript 21_gse141910_analysis.R
```

The workflow selects 166 DCM and 166 control samples, maps GeneIDs to gene symbols, retains the highest-mean-expression entry for duplicate symbols, applies expression filtering and TMM normalization, and fits a voom/limma group-only model with a DCM-minus-control contrast. It uses the full moderated-t ranking for enrichment analysis, with a fixed random seed of 42 and gene set size limits of 15–500.

Outputs are written to `GSE141910_output/`, including differential expression results, the full ranked list, Hallmark and KEGG enrichment results, and `sessionInfo.txt`.

## Reproducibility notes

The supplied ranked list was exported from preserved differential expression results. Preparing this repository did not involve rerunning differential expression or enrichment analysis. The supplied script documents the complete analysis workflow; an exact rerun of this packaged script has not been independently verified here.

Package and gene set database versions are not pinned in this repository. Changes in software versions or the gene sets retrieved by `msigdbr` may affect rerun results. The workflow records the R session information when executed.

SHA-256 checksums of the supplied analysis files:

```text
0f61821f00258a0623280e512164b0957801823905f7c13c2147273cbac34a92  21_gse141910_analysis.R
2c0eb5e85879218e8f1a5685a0b2d3f0a9e37807996202b085b87b1a17c32819  GSE141910_ranked_list.csv
```

## Data access

The clinical dataset is not deposited in this repository. Access to the clinical dataset is available from the corresponding author upon reasonable request, subject to approval by the Institutional Ethics Committee of West China Hospital. Public transcriptomic data are available under GEO accession GSE141910.
