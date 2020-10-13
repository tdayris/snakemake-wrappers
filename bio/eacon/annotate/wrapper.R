#!/bin/R

# This is the snakemake wrapper for EaCoN annotate

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

library(package = "EaCoN", quietly = TRUE);
library(package = "devtools", quietly = TRUE);

# Load user data
rds_path <- base::as.character(
  x = snakemake@input[["rds"]]
)

extra <- "";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["extra"]])
  );
}

grd_path <- base::normalizepath(
  path = base::as.character(x=snakemake@input[["grd"]])
);

base::Sys.setenv(
  PATH = paste(
    base::Sys.getenv("PATH"),
    grd_path,
    sep=":"
  )
);

# Build command line
command <- base::paste0(
  "EaCoN::ASCN.ff(",
  "RDS.file = rds_path, ",
  "grd = grd_path, "
  extra,
  ")"
);

# Run EaCoN
base::eval(
  base::parse(
    text = command
  )
);


base::Sys.setenv(
  PATH = paste(
    base::Sys.getenv("PATH"),
    snakemake@config[["params"]][["scripts"]],
    sep=":"
  )
);

EaCoN::Annotate.ff(
  RDS.file = snakemake@input[["rds"]],
  author.name = "STRonGR",
  ldb = snakemake@input[["ldb"]],
  solo = TRUE
);
