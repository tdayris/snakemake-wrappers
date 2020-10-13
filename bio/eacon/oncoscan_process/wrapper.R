#!/bin/R

# This is the snakemake wrapper for EaCoN OS.Process

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"


library(package = "EaCoN", quietly = TRUE);

atchc <- base::as.character(x = snakemake@input[["ATChannelCel"]]);
gcchc <- base::as.character(x = snakemake@input[["GCChannelCel"]]);

# Gather extra parameters
extra <- ", force=TRUE";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
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
  "EaCoN::OS.Process(",
  "ATChannelCel = atchc, ",
  "GCChannelCel = gcchc, ",
  "samplename = sample_name",
  extra,
  ")"
);
base::message(command);

base::eval(
  base::parse(
    text = command
  )
);
