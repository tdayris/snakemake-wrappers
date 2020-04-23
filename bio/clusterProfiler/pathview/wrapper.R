#!/usr/bin/R

# This script takes an gsea object from clusterProfiler
# and builds a GSEA-like plot

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Handle enrichment
base::library(package = "pathview", quietly = TRUE);
# Handle graphics
base::library(package = "Cairo", quietly = TRUE);

genelist <- base::readRDS(
  file = base::as.character(x = snakemake@input[["genelist"]])
);

pathway_id <- NULL;
if ("pathway_id" %in% base::names(snakemake@params)) {
  pathway_id <- base::as.character(snakemake@params[["pathway_id"]]);
} else {
  base::message("Missing parameter: pathway_id");
  quit(1);
}

extra <- "gene.data = genelist, pathway.id = pathway_id";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
}
base::message("Libraries and input data loaded");

command <- base::paste0(
  "pathview::pathview(",
  extra,
  ")"
);
base::message(command);

base::eval(
  base::parse(
    text = command
  )
);

base::file.rename(
  base::paste(pathway_id, "pathview", "png", sep="."),
  snakemake@output[["png"]]
);
