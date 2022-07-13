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
base::library(package = "org.Hs.eg.db", quietly = TRUE);
base::library(package = "org.Mm.eg.db", quietly = TRUE);

# Loading input datasets
geneList <- base::readRDS(file = snakemake@input[["rds"]]);
gmt_path <- base::as.character(x=snakemake@input[["gmt"]]);
gmt <- clusterProfiler::read.gmt(gmt_path);

org <- "org.Hs.eg.db";
if ("org" %in% base::names(snakemake@params)) {
  paramorg <- base::as.character(x=snakemake@params[["org"]]);
  if (paramorg == "Mm") {
    org <- "org.Mm.eg.db";
  }
}

base::message("Libraries and datasets loaded");

# Building command line
genes <- names(geneList);
universe <- NULL;
extra <- "gene=genes, TERM2GENE=gmt";

if ("universe" %in% base::names(snakemake@input)) {
  universe <- readRDS(file = snakemake@input[["universe"]]);
  extra <- base::paste(
    extra,
    "universe=universe",
    sep=", "
  );
}

if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep=", "
  );
}

command <- base::paste0(
  "clusterProfiler::GSEA(",
  extra,
  ")"
);
base::message(command);

# Performing DisGeNET enrichment
enriched_terms <- base::eval(
  base::parse(
    text = command
  )
);

# Saving results
base::saveRDS(
  obj = enriched_terms,
  file = snakemake@output[["rds"]]
);

utils::write.table(
  x = enriched_terms,
  file = snakemake@output[["tsv"]]
);


if ("org" %in% base::names(snakemake@params)) {
  readable_terms <- DOSE::setReadable(
    x=enriched_terms,
    OrgDb=base::as.character(snakemake@params[["org"]]),
    keyType="auto"
  );

  if ("readable_rds" %in% base::names(snakemake@output)) {
    base::saveRDS(
      obj = readable_terms,
      file = base::as.character(x=snakemake@output[["readable_rds"]])
    );
  }

  if ("readable_tsv" %in% base::names(snakemake@output)) {
    utils::write.table(
      x = readable_terms,
      file = base::as.character(x=snakemake@output[["readable_tsv"]])
    );
  }
}

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
