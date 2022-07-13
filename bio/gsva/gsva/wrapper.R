#!/usr/bin/R

# This script takes a GMT file and loads it into R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2021, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

base::library(package = "gsva", quietly=TRUE);

expr <- utils::read.table(
  file = base::as.character(x=snakemake@input[["expr"]]),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
);

gset <- base::readRDS(
  file = base::as.character(x=snakemake@input[["gset"]])
);
extra <- "expr = expr, gset.idx.list = gset";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
}
base::message("Libraries and input data loaded");

# Build command line
command <- base::paste0(
  "gsva::gsva(",
  heatplot_extra,
  ")"
);
base::message(command);



gsva_results <- base::eval(
  base::parse(
    text = command
  )
);

base::saveRDS(
  object = gsva_results,
  file = base::as.character(snakemake@output[["rds"]])
);

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
