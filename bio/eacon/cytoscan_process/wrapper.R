#!/bin/R

# This is the snakemake wrapper for EaCoN Process

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"


library(package = "EaCoN", quietly = TRUE);

cel <- base::as.character(x = snakemake@input[["cel"]]);

# Gather extra parameters
extra <- ", force=TRUE";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["extra"]])
  );
}

# Get sample name from wildcards
sample_name <- base::as.character(
  x = snakemake@wildcards["sample"]
);

# Build command line
command <- base::paste0(
  "EaCoN::CS.Process(CEL = cel, samplename = sample_name, ",
  extra,
  ")"
);
base::message(command);

base::eval(
  base::parse(
    text = command
  )
);
