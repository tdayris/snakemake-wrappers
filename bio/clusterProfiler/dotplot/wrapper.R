#!/usr/bin/R

# This script takes an enriched terms object from clusterProfiler
# and builds a dotplot of most enriched terms or provided list of
# pathways

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

# Handle enrichment
base::library(package = "clusterProfiler", quietly = TRUE);
# Handle graphics
base::library(package = "Cairo", quietly = TRUE);
base::library(package = "enrichplot", quietly = TRUE);

enriched <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
);

dotplot_extra <- "object = enriched";
if ("dotplot_extra" %in% base::names(snakemake@params)) {
  dotplot_extra <- base::paste(
    dotplot_extra,
    snakemake@params[["dotplot_extra"]],
    sep = ", "
  );
}
base::message("Libraries and input data loaded");

# Build command line
command <- base::paste0(
  "enrichplot::dotplot(",
  dotplot_extra,
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

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
