#!/bin/R

# This is the snakemake wrapper for EaCoN OS.Process

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");


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

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
