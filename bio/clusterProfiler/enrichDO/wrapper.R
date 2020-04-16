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

extra <- "gene = geneList, readable = FALSE, universe = names(geneList)";
if ("enrichDO_extra" %in% snakemake@params) {
  extra <- base::paste(
    extra,
    snakemake@params[["enrichDO_extra"]],
    sep = ", "
  )
}

command <- base::paste0(
  "clusterProfiler::enrichDO(",
  extra,
  ")"
);
base::message("Libraries and datasets loaded");
base::message(command);

# Performing DOSE enrichment
edo <- base::eval(
  base::parse(
    text = command
  )
);

edo <- clusterProfiler::setReadable(
    edo,
    organism,
    keytype = "ENTREZ"
);

# Saving results
base::saveRDS(
  obj = edo,
  file = snakemake@output[["rds"]]
);

utils::write.table(
  x = edo,
  file = snakemake@output[["tsv"]]
);
