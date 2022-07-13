#!/usr/bin/R

# This script takes a IsoformSwitchAnalyseR object
# and tests exons with differential expression

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# The main package
base::library(package = "IsoformSwitchAnalyzeR", quietly = TRUE);
base::library(package = "DEXSeq", quietly = TRUE);

# Gathering input dataset and parameters
switch_list <- base::readRDS(
  file = base::as.character(x = snakemake@input[["switch_list"]])
);

extra <- "switchAnalyzeRlist = switch_list";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
      extra,
      snakemake@params[["extra"]],
      sep = ", "
  );
}

command <- base::paste0(
  "IsoformSwitchAnalyzeR::isoformSwitchTestDEXSeq(",
  extra,
  ")"
);
base::message("Libraries and input data loaded");
base::message(command);

# Running command
tested <- base::eval(
  base::parse(
    text = command
  )
);

# Saving results
if ("rds" %in% base::names(snakemake@output)) {
  base::saveRDS(
    obj = tested,
    file = snakemake@output[["rds"]]
  );
}

if ("tsv" %in% base::names(snakemake@output)) {
  tested_table <- base::as.data.frame(
    x = tested$isoformSwitchAnalysis,
    stringsAsFactors = FALSE
  );

  utils::write.table(
    x = tested_table,
    file = base::as.character(snakemake@output[["tsv"]]),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  );
}
