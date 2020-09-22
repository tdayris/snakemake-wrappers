#!/usr/bin/R

# This script takes a IsoformSwitchAnalyseR object
# and extract sequences from targets

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# The main package
base::library(package = "IsoformSwitchAnalyzeR", quietly = TRUE);

# Gathering input dataset and parameters
switch_list <- base::readRDS(
  file = base::as.character(x = snakemake@input[["switch_list"]])
);

prefix <- "isoformSwitchAnalyzeR_isoform";
if ("prefix" %in% base::names(snakemake@params)) {
  prefix <- base::as.character(snakemake@params[["prefix"]]);
} else if ("fasta" %in% base::names(snakemake@output)) {
  prefix <- base::gsub("_nt.fasta", snakemake@output[["fasta"]])
} else if ("aa_sequence" %in% base::names(snakemake@output)) {
  prefix <- base::gsub("_AA.fasta", snakemake@output[["fasta"]])
}

extra <- "switchAnalyzeRlist = switch_list, outputPrefix = prefix";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
      extra,
      snakemake@params[["extra"]],
      sep = ", "
  );
}

command <- base::paste0(
  "IsoformSwitchAnalyzeR::extractSequence(",
  extra,
  ")"
);
base::message("Libraries and input data loaded");
base::message(command);

# Running command
extracted <- base::eval(
  base::parse(
    text = command
  )
);

# Saving results
base::saveRDS(
  obj = extracted,
  file = snakemake@output[["rds"]]
);
