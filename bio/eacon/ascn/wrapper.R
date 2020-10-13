#!/bin/R

# This is the snakemake wrapper for EaCoN ASCN

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

library(package = "EaCoN", quietly = TRUE);

# Load user data
rds_path <- base::as.character(
  x = snakemake@input[["rds"]]
)

extra <- "";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste0(
    ",",
    base::as.character(x = snakemake@params[["extra"]])
  );
}

# Build command line
command <- base::paste0(
  "EaCoN::ASCN.ff(",
  "RDS.file = rds_path",
  extra,
  ")"
);

# Run EaCoN
base::eval(
  base::parse(
    text = command
  )
);
