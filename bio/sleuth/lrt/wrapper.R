#!/bin/R

# Snakemake wrapper for sleuth lrt
# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

base::library(package = "sleuth", quietly = TRUE);

threads <- snakemake@threads;
mc.cores <- threads;
options(mc.cores = threads);

so <- base::readRDS(
    file = snakemake@input[["rds"]]
);

extra <- "";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::as.character(
    x = snakemake@params[["extra"]]
  );
}

command <- base::paste(
  "sleuth::sleuth_lrt(",
  "obj = so",
  extra,
  ")",
  sep = ", "
);

## Finally run sleuth prep
so_lrt <- base::eval(
  base::parse(
    text = command
  )
);


base::saveRDS(
  object = so_lrt,
  file = snakemake@output[["rds"]]
);
