#!/usr/bin/R

# This script takes a geneList object and performs
# an enrichment analysis based on DOSE database

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

extra <- "geneList, OrgDb = organism";
if ("gseGO" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["gseGO"]],
    sep = ", "
  )
}

command <- base::paste0(
  "clusterProfiler::gseGO(",
  extra,
  ")"
);
base::message("Libraries and datasets loaded");
base::message(command);

# Performing DOSE enrichment
gse_go <- base::eval(
  base::parse(
    text = command
  )
);
base::message("GSEA done");

gse_go <- clusterProfiler::setReadable(
    x = gse_go,
    OrgDb = organism
);

# Saving results
base::saveRDS(
  obj = gse_go,
  file = snakemake@output[["rds"]]
);

utils::write.table(
  x = gse_go,
  file = snakemake@output[["tsv"]]
);
