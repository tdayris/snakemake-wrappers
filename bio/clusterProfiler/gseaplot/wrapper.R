#!/usr/bin/R

# This script takes an gsea object from clusterProfiler
# and builds a GSEA-like plot

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Handle enrichment
base::library(package = "clusterProfiler", quietly = TRUE);
# Handle graphics
base::library(package = "Cairo", quietly = TRUE);
base::library(package = "enrichplot", quietly = TRUE);

gsea <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
);

extra <- "x = gsea";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
} else {
  extra <- base::paste(
    extra,
    "geneSetID = 1",
    sep = ", "
  )
}
base::message("Libraries and input data loaded");

command <- base::paste0(
  "enrichplot::gseaplot2(",
  extra,
  ")"
);
base::message(command);

# Build plot
png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px",
  type = "cairo"
);

base::eval(
  base::parse(
    text = command
  )
);

dev.off()
