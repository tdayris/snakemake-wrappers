#!/usr/bin/R

# This script takes a geneList object and performs
# an classification based on GO

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Perform gene enrichment
base::library(package = "clusterProfiler", quietly = TRUE);
# Loading databases
base::library(package = "org.Hs.eg.db", quietly = TRUE);
base::library(package = "org.Mm.eg.db", quietly = TRUE);


# Loading input dataset
geneList <- base::readRDS(
  file = snakemake@input[["rds"]]
);

organism <- org.Hs.eg.db;
if ("organism" %in% base::names(snakemake@params)) {
  if (snakemake@params[["organism"]] == "Mm") {
    organism <- org.Mm.eg.db;
  }
}

extra <- "gene = base::names(geneList), OrgDb = organism, readable = TRUE";
if ("groupGO_extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["groupGO_extra"]],
    sep = ", "
  )
}

command <- base::paste0(
  "clusterProfiler::groupGO(",
  extra,
  ")"
);
base::message("Libraries and datasets loaded");
base::message(command);

# Performing DOSE enrichment
gse_do <- base::eval(
  base::parse(
    text = command
  )
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
