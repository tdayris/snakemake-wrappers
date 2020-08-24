#!/usr/bin/R

# This script takes an list of Salmon directories
# and builds a IsoformSwitchAnalyseR object

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"


# The main package
base::library(package = "IsoformSwitchAnalyzeR", quietly = TRUE);

# Gathering input files
quant_files <- sapply(
  snakemake@input[["quant"]],
  function(quant) base::as.character(x = quant)
);

# Guessing the argument type in importIsoformExpression
import_type <- ifelse(
  all(file_test("-f", quant_files)),
  "sampleVector",
  "parentDir"
);


# Building command line arguments
extra <- base::paste0(import_type, " = quant_files");
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
}

# Buiding command line itself
command <- base::paste0(
  "IsoformSwitchAnalyzeR::importIsoformExpression(",
  extra,
  ")"
);
base::message("Libraries and input data loaded");
base::message(command);

# Running command
isar <- base::eval(
  base::parse(
    text = command
  )
);

# Saving results
base::saveRDS(
  obj = isar,
  file = snakemake@output[["rds"]]
);
