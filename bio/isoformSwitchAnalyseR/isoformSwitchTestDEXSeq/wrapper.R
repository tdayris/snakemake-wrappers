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
base::saveRDS(
  obj = tested,
  file = snakemake@output[["rds"]]
);
