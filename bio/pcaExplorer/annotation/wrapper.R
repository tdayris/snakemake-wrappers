#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script takes a tsv formatted text file
# composed of the following columns:
# gene_id, transcript_id, gene_name
# And builds an object designed for pcaExplorer

# Load libraries
base::library(package = "DelayedArray");
base::library(package = "IRanges");
base::library(package = "readr");

# Load datasets
dds_path <- base::as.character(
  x = snakemake@input[["dds"]]
);
dds <- base::readRDS(file = dds_path);

tx2gene_path <- base::as.character(
  x = snakemake@input[["tr2gene"]]
);
tx2gene <- utils::read.table(
  file = tx2gene_path,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
);

# Remove un-used transcript_id
tx2gene <- tx2gene[, c("V1", "V3")];
IRanges::colnames(tx2gene) <- c("gene_id", "gene_name");
tx2gene <- DelayedArray::unique(tx2gene);
base::row.names(tx2gene) <- tx2gene$gene_id;

# Build final dataframe
gene_names <- tx2gene[base::row.names(x = dds), ];
gene_names$gene_id <- NULL;
annotation <- base::data.frame(
  gene_name = gene_names,
  row.name = IRanges::rownames(x = dds),
  stringsAsFactors = FALSE
);

annot_output <- base::as.character(
  x = snakemake@output[["annotation"]]
);
base::saveRDS(
  object = annotation,
  file = annot_output
);
