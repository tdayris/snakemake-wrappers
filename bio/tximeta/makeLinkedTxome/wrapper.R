#!/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Loading library
base::library(package = "tximeta", quietly = TRUE);
base::message("Libraries loaded");

# Defining output directory
bioc_dir <- base::as.character(x = snakemake@output[["bioc_dir"]]);
base::dir.create(bioc_dir, showWarnings = TRUE);
tximeta::setTximetaBFC(dir = bioc_dir, quiet = FALSE);
write("Cache destination set");

# Loading user defined information about the genome
index <- base::as.character(x = snakemake@input[["index"]]);
source <- base::as.character(x = snakemake@params[["source"]]);
organism <- base::as.character(x = snakemake@params[["organism"]]);
release <- base::as.character(x = snakemake@params[["release"]]);
genome <- base::as.character(x = snakemake@params[["genome"]]);

gtf <- NULL;
if ("gtf" %in% base::names(snakemake@input)) {
  gtf <- base::as.character(x = snakemake@input[["gtf"]]);
} else if ("gtf" %in% base::names(snakemake@params)) {
  gtf <- base::as.character(x = snakemake@params[["gtf"]]);
}

fasta <- NULL;
if ("fasta" %in% base::names(snakemake@input)) {
  fasta <- base::as.character(x = snakemake@input[["fasta"]]);
} else if ("fasta" %in% base::names(snakemake@params)) {
  fasta <- base::as.character(x = snakemake@params[["fasta"]]);
}

write <- FALSE;
json_path <- NULL;
if ("json" %in% base::names(snakemake@output)) {
  json_path <- base::as.character(x = snakemake@output[["json"]]);
  write <- TRUE;
}
base::message("Parametres set");

tximeta::makeLinkedTxome(
  indexDir = index,
  source = source,
  organism = organism,
  release = release,
  genome = genome,
  fasta = fasta,
  gtf = gtf,
  write = write,
  json = json_path
);
base::message("Process over");
