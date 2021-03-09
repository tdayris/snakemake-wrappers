#!/usr/bin/R

# This script takes a geneList object and performs
# an enrichment analysis based on GO database

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

# Perform gene enrichment
base::library(package = "DOSE", quietly = TRUE);
base::library(package = "clusterProfiler", quietly = TRUE);
# Loading databases
base::library(package = "org.Hs.eg.db", quietly = TRUE);

# Loading input dataset
geneList <- base::readRDS(
  file = snakemake@input[["rds"]]
);

extra <- "geneList, verbose = TRUE";
if ("gseDO_extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["gseDO_extra"]],
    sep = ", "
  )
}

command <- base::paste0(
  "DOSE::gseDO(",
  extra,
  ")"
);
base::message("Libraries and datasets loaded");
base::message(command);

# Performing GO enrichment
gse_do <- base::eval(
  base::parse(
    text = command
  )
);

gse_do <- DOSE::setReadable(
    x = gse_do,
    OrgDb = org.Hs.eg.db
);

# Saving results
base::saveRDS(
  obj = gse_do,
  file = snakemake@output[["rds"]]
);

utils::write.table(
  x = gse_do,
  file = snakemake@output[["tsv"]]
);

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
