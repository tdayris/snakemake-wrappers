#!/bin/R

# This is the snakemake wrapper for EaCoN segmentation

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

base::library(package = "EaCoN", quietly = TRUE);

input_rds <- base::as.character(
  x = snakemake@input[["rds"]]
);

extra <- "force=TRUE";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::as.character(x = snakemake@params[["extra"]]);
}

command <- base::paste0(
  "EaCoN::Segment.ff(RDS.file = input_rds, ",
  extra,
  ")"
);
base::message(command);

base::eval(
  base::parse(
    text=command
  )
);

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
