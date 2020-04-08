#!/usr/bin/R

# This script takes a deseq2 dataset object and estimates
# size factors for further normalization

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Differential Gene expression
base::library(package = "DESeq2", quietly = TRUE);

# Cast input path as character
dds_path <- base::as.character(x = snakemake@input[["dds"]]);
dds <- base::readRDS(dds_path);


# Check if user provided optional parameters
extra <- "";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["extra"]])
  );
}


# Create object
dds <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::estimateSizeFactors(dds", extra, ");"
    )
  )
);

# Save as RDS
output_path <- base::as.character(snakemake@output[["dds"]]);
base::saveRDS(
  obj = dds,
  file = output_path
);
