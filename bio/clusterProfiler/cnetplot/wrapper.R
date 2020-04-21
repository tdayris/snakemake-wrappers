#!/usr/bin/R

# This script takes an enriched terms object from clusterProfiler
# and builds a cnetplot of most enriched terms or provided list of
# pathways

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Handle enrichment
base::library(package = "clusterProfiler", quietly = TRUE);
# Handle graphics
base::library(package = "Cairo", quietly = TRUE);
base::library(package = "enrichplot", quietly = TRUE);

enriched <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
);

fc <- NULL;
if ("gene_list" %in% base::names(snakemake@input)) {
  fc <- readRDS(
    file = base::as.character(snakemake@input[["gene_list"]])
  );
}

cnetplot_extra <- "x = enriched, foldChange = fc";
if ("cnetplot_extra" %in% base::names(snakemake@params)) {
  cnetplot_extra <- base::paste(
    cnetplot_extra,
    snakemake@params[["cnetplot_extra"]],
    sep = ", "
  );
}
base::message("Libraries and input data loaded");

# Build command line
command <- base::paste0(
  "enrichplot::cnetplot(",
  cnetplot_extra,
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
