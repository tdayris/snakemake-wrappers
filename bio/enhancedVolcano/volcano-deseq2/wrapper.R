#!/usr/bin/R

# This script builds a volcano plot from a deseq2 tsv file

# Loading libraries
base::library(package = "readr", quietly = TRUE);  # Read large dataset
base::library(package = "EnhancedVolcano", quietly = TRUE);  # Perform volcano plot

dataset <- utils::read.table(
  file = snakemake@input[["deseq2_tsv"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
);


padj_threshold <- 0.05
if ("padj_threshold" %in% snakemake@params) {
  padj_threshold <- snakemake@params[["padj_threshold"]];
}

fc_threshold <- 0.5
if ("fc_threshold" %in% snakemake@params) {
  fc_threshold <- snakemake@params[["fc_threshold"]];
}


grDevices::png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px"
);

EnhancedVolcano::EnhancedVolcano(
  dataset,
  lab = rownames(dataset),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = padj_threshold,
  FCcutoff = fc_threshold,
  title = "Volcano plot",
  subtitle = snakemake@params[["condition_name"]],
  pLabellingCutoff = 10^-1000000
);


dev.off();
