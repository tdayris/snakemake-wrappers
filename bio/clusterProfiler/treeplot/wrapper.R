#!/usr/bin/R

# This script takes an enriched terms object from clusterProfiler
# and builds a tree of most enriched terms or provided list of
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
# base::library(package = "Cairo", quietly = TRUE);
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

extra_paiwise_termism <- "x = enriched";
if ("extra_paiwise_termism" %in% base::names(snakemake@params)) {
  extra_paiwise_termism <- base::paste(
    extra_paiwise_termism,
    snakemake@params[["extra_paiwise_termism"]],
    sep = ", "
  );
}

extra_treeplot <- "x = enriched";
if ("extra_treeplot" %in% base::names(snakemake@params)) {
  extra_treeplot <- base::paste(
    extra_treeplot,
    snakemake@params[["extra_treeplot"]],
    sep = ", "
  );
}
base::message("Libraries and input data loaded");

# Build paiwise_termism command line
command <- base::paste0(
  "enrichplot::pairwise_termsim(",
  extra_paiwise_termism,
  ")"
);
base::message(command);

# Execute pairwise_termsim clustering
pairs <- base::eval(base::parse(text = command));

# Build treeplot command line
command <- base::paste0(
  "enrichplot::treeplot(",
  extra_treeplot,
  ")"
);
base::message(command);

# Build plot
png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px",
  # type = "cairo"
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
