#!/usr/bin/env Rscript

# GSE141910 transcriptomic analysis
#
# Inputs expected in ./GSE141910_input:
# - GSE141910_raw_counts_GRCh38.p13_NCBI.tsv.gz
# - Human.GRCh38.p13.annot.tsv.gz
#
# Outputs are written to ./GSE141910_output. The script retrieves the public
# GSE141910 phenotype table from GEO but does not download expression files.

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(dplyr)
})

set.seed(42)

input_dir <- "GSE141910_input"
output_dir <- "GSE141910_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(
  input_dir, "GSE141910_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)
annotation_path <- file.path(input_dir, "Human.GRCh38.p13.annot.tsv.gz")
stopifnot(file.exists(counts_path), file.exists(annotation_path))

gse141910 <- getGEO("GSE141910", GSEMatrix = TRUE, getGPL = FALSE)
pheno141910 <- pData(gse141910[[1]])

raw_counts <- read.delim(
  counts_path,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
annotation <- read.delim(
  annotation_path,
  header = TRUE,
  check.names = FALSE
)

pheno141910$group <- case_when(
  grepl("Dilated cardiomyopathy|DCM", pheno141910$characteristics_ch1) ~ "DCM",
  grepl("Non-Failing|Donor", pheno141910$characteristics_ch1) ~ "Control",
  TRUE ~ "Other"
)
pheno_use <- pheno141910[
  pheno141910$group %in% c("DCM", "Control"),
  ,
  drop = FALSE
]
common <- intersect(colnames(raw_counts), pheno_use$geo_accession)
counts_use <- raw_counts[, common, drop = FALSE]
pheno_use <- pheno_use[match(common, pheno_use$geo_accession), , drop = FALSE]
stopifnot(
  length(common) == 332L,
  sum(pheno_use$group == "DCM") == 166L,
  sum(pheno_use$group == "Control") == 166L
)

gene_map <- setNames(annotation$Symbol, as.character(annotation$GeneID))
mapped_symbols <- gene_map[rownames(counts_use)]
valid <- !is.na(mapped_symbols) & mapped_symbols != "" & mapped_symbols != "-"
counts_use <- counts_use[valid, , drop = FALSE]
mapped_symbols <- mapped_symbols[valid]

counts_use <- as.matrix(counts_use)
mode(counts_use) <- "numeric"
deduplication_map <- data.frame(
  symbol = mapped_symbols,
  mean_expression = rowMeans(counts_use),
  row_index = seq_along(mapped_symbols),
  stringsAsFactors = FALSE
) %>%
  group_by(symbol) %>%
  slice_max(mean_expression, n = 1, with_ties = FALSE) %>%
  ungroup()
counts_use <- counts_use[deduplication_map$row_index, , drop = FALSE]
rownames(counts_use) <- deduplication_map$symbol

dge <- DGEList(counts = counts_use, group = pheno_use$group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~0 + factor(pheno_use$group))
colnames(design) <- c("Control", "DCM")
voom_fit <- voom(dge, design, plot = FALSE)
contrast <- makeContrasts(DCM_vs_Control = DCM - Control, levels = design)
fit <- lmFit(voom_fit, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

de_results <- topTable(fit, coef = 1, number = Inf, sort.by = "none")
de_results$gene_symbol <- rownames(de_results)
write.csv(
  de_results,
  file.path(output_dir, "DE_results_GSE141910.csv"),
  row.names = FALSE
)

gene_ranks <- setNames(de_results$t, de_results$gene_symbol)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)
ranked_list <- data.frame(
  gene_symbol = names(gene_ranks),
  t = as.numeric(gene_ranks),
  stringsAsFactors = FALSE
)
write.csv(
  ranked_list,
  file.path(output_dir, "GSE141910_ranked_list.csv"),
  row.names = FALSE
)

hallmark_membership <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_pathways <- split(
  hallmark_membership$gene_symbol,
  hallmark_membership$gs_name
)

kegg_membership <- tryCatch(
  msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:KEGG"
  ),
  error = function(e) {
    c2_membership <- msigdbr(species = "Homo sapiens", collection = "C2")
    kegg_subcollections <- unique(c2_membership$gs_subcollection)[
      grepl("KEGG", unique(c2_membership$gs_subcollection))
    ]
    if (length(kegg_subcollections)) {
      msigdbr(
        species = "Homo sapiens",
        collection = "C2",
        subcollection = kegg_subcollections[[1]]
      )
    } else {
      c2_membership[grepl("^KEGG_", c2_membership$gs_name), , drop = FALSE]
    }
  }
)
kegg_pathways <- split(kegg_membership$gene_symbol, kegg_membership$gs_name)

hallmark_results <- fgsea(
  pathways = hallmark_pathways,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  nPermSimple = 10000
)
kegg_results <- fgsea(
  pathways = kegg_pathways,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  nPermSimple = 10000
)

prepare_fgsea_export <- function(results) {
  output <- as.data.frame(results)
  output$leadingEdge <- vapply(
    output$leadingEdge,
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
  )
  output[order(output$padj, output$pval, output$pathway), , drop = FALSE]
}

write.csv(
  prepare_fgsea_export(hallmark_results),
  file.path(output_dir, "GSEA_Hallmark_all.csv"),
  row.names = FALSE
)
write.csv(
  prepare_fgsea_export(kegg_results),
  file.path(output_dir, "GSEA_KEGG_all.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)
