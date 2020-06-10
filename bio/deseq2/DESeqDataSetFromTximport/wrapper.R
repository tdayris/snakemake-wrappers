#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script takes a tximport object and builds a deseq2 dataset
# for each formula given to snakemake.

# Perform actual count importation in R
base::library(package = "tximport", quietly = TRUE);
# Read faster!
base::library(package = "readr", quietly = TRUE);
# Importing inferential replicates
base::library(package = "jsonlite", quietly = TRUE);
# Differential Gene expression
base::library(package = "DESeq2", quietly = TRUE);

base::write("Libraries loaded.", stderr());

# Load txi object
txi_rds_path <- base::as.character(x = snakemake@input[["tximport"]]);
txi <- base::readRDS(
  file = txi_rds_path
);

# Load experimental design
coldata_path <- base::as.character(x = snakemake@input[["coldata"]]);
coldata <- utils::read.table(
  file = coldata_path,
  sep = "\t",
  header = TRUE
);
rownames(coldata) <- coldata$Sample_id;

count_filter <- 0;
if ("count_filter" %in% names(snakemake@params)) {
  count_filter <- base::as.numeric(
    x = snakemake@params[["count_filter"]]
  );
}
keep <- rowSums(counts(dds)) > count_filter;
dds <- dds[keep, ];

# Cast formula as formula instead of string
formula <- stats::as.formula(
  object = snakemake@params[["design"]]
);

base::message("Input dataset and options recovered.");

# Create dds object
dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi,
  colData = coldata,
  design = formula
);
base::write("DESeqDataSet built.", stderr());



# Save as RDS
output_path <- base::as.character(x = snakemake@output[["dds"]]);
base::saveRDS(
  obj = dds,
  file = output_path
);

base::write("Process over.", stderr());
