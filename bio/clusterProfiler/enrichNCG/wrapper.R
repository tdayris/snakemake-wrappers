#!/usr/bin/R

# This script takes a geneList object and performs
# an enrichment analysis based on NCG database

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Perform gene enrichment
base::library(package = "DOSE", quietly = TRUE);
base::library(package = "clusterProfiler", quietly = TRUE);
# Loading databases
base::library(package = "org.Hs.eg.db", quietly = TRUE);


# Loading input dataset
geneList <- base::readRDS(
  file = snakemake@input[["rds"]]
);

extra <- "gene = base::names(geneList), readable = TRUE";
if ("enrichNCG_extra" %in% snakemake@params) {
  extra <- base::paste(
    extra,
    snakemake@params[["enrichNCG_extra"]],
    sep = ", "
  )
}

command <- base::paste0(
  "DOSE::enrichNCG(",
  extra,
  ")"
);
base::message("Libraries and datasets loaded");
base::message(command);

# Performing DOSE enrichment
encg <- base::eval(
  base::parse(
    text = command
  )
);

# Saving results
base::saveRDS(
  obj = encg,
  file = snakemake@output[["rds"]]
);

utils::write.table(
  x = encg,
  file = snakemake@output[["tsv"]]
);
