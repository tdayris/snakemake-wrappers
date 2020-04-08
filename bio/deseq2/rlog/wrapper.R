#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script takes a deseq2 dataset object and performs
# a rlog transformation on it

base::library("DESeq2");                 # Differential Gene expression
base::library("SummarizedExperiment");   # Handle large datasets

# Cast input path as character
dds_path <- base::as.character(x = snakemake@input[["dds"]]);
dds <- base::readRDS(file = dds_path);

# Recover extra parameters
extra <- "";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
    "object = dds",
    base::as.character(x = snakemake@params[["extra"]]),
    sep = ", "
  );
}

# Create object
rlog <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::rlogTransformation(", extra, ");"
    )
  )
);

tsv <- SummarizedExperiment::assay(rlog);

# Save results
output_rds <- base::as.character(snakemake@output[["rds"]]);
base::saveRDS(
  obj = rlog,
  file = output_rds
);


output_table <- base::as.character(snakemake@output[["tsv"]]);
utils::write.table(
  x = tsv,
  file = output_table,
  sep = "\t",
  quote = FALSE
);
