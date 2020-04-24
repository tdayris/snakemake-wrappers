#!/usr/bin/R

# This script takes an list of Salmon directories
# and builds a IsoformSwitchAnalyseR object

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"


# The main package
base::library(package = "isoformSwitchAnalyzeR", quietly = TRUE);

# Gathering parameters
quant_files <- sapply(
  snakemake@input[["quant"]],
  function(quant) base::as.character(x = quant)
);


extra <- "parentDir = quant_files";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
}

# Buiding command line itself
command <- base::paste0(
  "InsoformSwitchAnalyseR::importIsoformExpression(",
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
