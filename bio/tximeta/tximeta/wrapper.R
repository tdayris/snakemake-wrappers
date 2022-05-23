#!/usr/bin/env R

# Snakemake wrapper for tximeta main function
# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2022, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file <- base::file(snakemake@log[[1]], open="wt")
base::sink(file=log.file);
base::sink(file=log.file, type="message");

base::library(packages="tximeta");


if ("json" %in% base::names(snakemake@input)) {
  tximeta::loadLinkedTxome(
    jsonFile=base::as.character(x=snakemake@input[["json"]])
  );
}

# Get user defined optional arguments
parameters <- 'coldata=snakemake@input[["coldata"]]';
if ("extra" %in% base::names(snakemake@params)) {
  parameters <- base::paste(
    parameters,
    base::as.character(x=snakemake@params[["extra"]]),
    sep=", "
  );
}

# Build and execute command line with user defined arguments
tximeta_cmd <- base::paste0("tximeta::tximeta(", parameters, ")");
base::message("tximeta command line:", tximeta_cmd);
se <- base::eval(base::paste(text=tximeta_cmd));

# Save raw RDS on user request
if ("rds" %in% base::names(snakemake@output)) {
  base::saveRDS(
    object=se,
    file=base::as.character(x=snakemake@output[["rds"]])
  );
}

# Compute and save counts-to-gene aggregation on user request
if "gene_rds" %in% base::names(snakemake@output) {

  # Build command line, with user defined optional arguments
  parameters <- 'object=se';
  if ("summarize_to_gene_extra" %in% base::names(snakemake@params)) {
    parameters <- base::paste(
      parameters,
      base::as.character(x=snakemake@params[["extra"]]),
      sep=", "
    );
  }
  summarize_cmd <- base::paste0("tximeta::summarizeToGene(", parameters, ")");
  base::message("summarizeToGene command line:", summarize_cmd);

  # Execute and save results
  gene_se <- base::eval(base::paste(text=tximeta_cmd));
  base::saveRDS(
    object=gene_se,
    file=base::as.character(x=snakemake@output[["gene_rds"]])
  );
}

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
