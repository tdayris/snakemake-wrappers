#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for maftools read.maf

# Many libraries are useless. This wrapper will be splitted in at least four
# of them.
base::library(package = "maftools", quietly = TRUE);

# Build command line:
maf <- base::as.character(x = snakemake@input[["maf"]]);

cmd_args <- "maf = maf";
if ("clinical_data" %in% base::names(snakemake@input)) {
  clinical_data_arg <- base::paste(
    "clinicalData",
    "=",
    base::as.character(x = snakemake@input[["clinical_data"]])
  );
  cmd_args <- base::paste(cmd_args, clinical_data_arg, sep=", ");
}


if ("extra" %in% base::names(snakemake@params)) {
  cmd_args <- base::paste(cmd_args, snakemake@params[["extra"]], sep=", ");
}

cmd_line <- base::paste0("maftools::read.maf(", cmd_args, ")");

# Executing maftools command
base::message(cmd_line);
maf_obj <- base::eval(base::parse(text = cmd_line));


# Saving results on demand
if ("rds" %in% base::names(snakemake@output)) {
  base::saveRDS(
    file = base::as.character(x = snakemake@output[["rds"]]),
    obj = maf_obj
  );
}

if ("summary_prefix" %in% base::names(snakemake@params)) {
  maftools:::write.mafSummary(
    maf = maf_obj,
    basename = base::as.character(x = snakemake@params[["summary_prefix"]])
  );
}
