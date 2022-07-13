#!/bin/R

# This is the snakemake wrapper for EaCoN annotate

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

library(package = "EaCoN", quietly = TRUE);
library(package = "devtools", quietly = TRUE);

# Load user data
rds_path <- base::as.character(
  x = snakemake@input[["rds"]]
)

grd_path <- base::normalizePath(
  path = base::as.character(x=snakemake@input[["grd"]])
);

ldb <- base::normalizePath(
  path = base::as.character(x=snakemake@input[["ldb"]])
);

base::Sys.setenv(
  PATH = paste(
    base::Sys.getenv("PATH"),
    base::dirname(grd_path),
    sep=":"
  )
);

print(base::Sys.getenv("PATH"));


EaCoN::Annotate.ff(
  RDS.file = rds_path,
  author.name = "BiGR",
  ldb = ldb,
  solo = TRUE
);


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
