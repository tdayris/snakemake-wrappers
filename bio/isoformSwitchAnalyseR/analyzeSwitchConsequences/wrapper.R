#!/usr/bin/R

# This script takes a IsoformSwitchAnalyseR object
# and add CPAT results

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

consequences <- sapply(
  snakemake@params[["consequences"]],
  function(consequence) base::as.character(x = consequence)
);

extra <- base::paste(
  "switchAnalyzeRlist = switch_list",
  "consequencesToAnalyze = consequences",
  sep = ", "
);
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
}

# Buiding command line itself
command <- base::paste0(
  "InsoformSwitchAnalyseR::analyzeSwitchConsequences(",
  extra,
  ")"
);
base::message("Libraries and input data loaded");
base::message(command);

# Running command
consequences <- base::eval(
  base::parse(
    text = command
  )
);

# Saving results
base::saveRDS(
  obj = consequences,
  file = snakemake@output[["rds"]]
);
