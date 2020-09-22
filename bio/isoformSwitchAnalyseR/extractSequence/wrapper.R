#!/usr/bin/R

# This script takes a IsoformSwitchAnalyseR object
# and extract sequences from targets

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# The main package
base::library(package = "isoformSwitchAnalyzeR", quietly = TRUE);

# Gathering input dataset and parameters
switch_list <- base::readRDS(
  file = base::as.character(x = snakemake@input[["switch_list"]])
);

extra <- "switchAnalyzeRlist = switch_list";
if ("extra" %in% base::anems(snakemake@params)) {
  extra <- base::paste(
      extra,
      snakemake@params[["extra"]],
      sep = ", "
  );
}

command <- base::paste0(
  "InsoformSwitchAnalyseR::extractSequence(",
  extra,
  ")"
);
base::message("Libraries and input data loaded");
base::message(command);

# Running command
base::eval(
  base::parse(
    text = command
  )
);

# Saving results
base::saveRDS(
  obj = edgn,
  file = snakemake@output[["rds"]]
);
