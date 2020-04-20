#!/usr/bin/R

# This script takes a geneList object and performs
# a GSEA on DisGeNET

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Perform gene enrichment
base::library(package = "clusterProfiler", quietly = TRUE);
base::library(package = "DOSE", quietly = TRUE);
# Loading databases
base::library(package = "org.Hs.eg.db", quietly = TRUE);


# Loading input dataset
geneList <- base::readRDS(
  file = snakemake@input[["rds"]]
);

extra <- "geneList, verbose = TRUE";
if ("gseDGN_extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["gseDGN_extra"]],
    sep = ", "
  )
}

command <- base::paste0(
  "DOSE::gseDGN(",
  extra,
  ")"
);
base::message("Libraries and datasets loaded");
base::message(command);

# Performing DOSE enrichment
gse_dgn <- base::eval(
  base::parse(
    text = command
  )
);

gse_dgn <- DOSE::setReadable(
    x = gse_dgn,
    OrgDb = org.Hs.eg.db
);

# Saving results
base::saveRDS(
  obj = gse_dgn,
  file = snakemake@output[["rds"]]
);

utils::write.table(
  x = gse_dgn,
  file = snakemake@output[["tsv"]]
);
