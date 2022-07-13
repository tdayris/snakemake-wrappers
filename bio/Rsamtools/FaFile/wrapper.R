#!/usr/bin/R

# Load a fasta file and return a RDS object
library(package="Rsamtools", quietly=TRUE);

# Load dataset
fasta <- Rsamtools::FaFile(
  file = base::as.character(x = snakemake@input[["fasta"]]),
  index = base::as.character(x = snakemake@input[["fai"]])
);

# Save as RDS object
base::saveRDS(
  object = fasta,
  file = base::as.character(x = snakemake@output[["rds"]])
);
