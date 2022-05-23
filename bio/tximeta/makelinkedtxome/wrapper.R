#!/use/bin/env R

# Snakemake wrapper for tximeta::makeLinkedTxome
# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2022, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
# log_file <- base::file(snakemake@log[[1]], open="wt");
# base::sink(file=log_file);
# base::sink(file=log_file, type="message");

base::library(package="tximeta", quietly=TRUE);

tximeta::makeLinkedTxome(
  indexDir=base::as.character(x=snakemake@input[["index"]]),
  source=base::as.character(x=snakemake@params[["source"]]),
  organism=base::as.character(x=snakemake@params[["organism"]]),
  release=base::as.character(x=snakemake@params[["release"]]),
  genome=base::as.character(x=snakemake@params[["genome"]]),
  fasta=base::as.character(x=snakemake@input[["fasta"]]),
  gtf=base::as.character(x=snakemake@input[["gtf"]]),
  write=TRUE,
  jsonFile=base::as.character(x=snakemake@output[[1]])
);


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
# base::sink(type="message");
# base::sink();
