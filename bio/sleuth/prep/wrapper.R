#!/bin/R

# Snakemake wrapper for sleuth prep
# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

base::library(package = "sleuth", quietly = TRUE);

threads <- snakemake@threads;
mc.cores <- threads;
options(mc.cores = threads);

command <- "sleuth::sleuth_prep"

# Load transcript to gene translation table if provided
target_mapping <- NULL;
if ("target_mapping" %in% names(snakemake@input)) {
  target_mapping <- utils::read.table(
    file = snakemake@input[["target_mapping"]],
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = "\t"
  );
  command <- base::paste(
    command,
    "target_mapping = target_mapping",
    sep = ", "
  );
}


# Load experimental design to perform differential analysis
meta <- utils::read.table(
    file = snakemake@input[["meta"]],
    header = TRUE,
    sep = "\t"
);
command <- base::paste(
  command,
  "sample_to_covariates = meta",
  sep = ", "
);


# Load user defined arguments
if ("extra" %in% names(snakemake@params)) {
  if (snakemake@params[["extra"]] != "") {
    command <- base::paste(
      command,
      snakemake@params[["extra"]],
      ")"
      sep = ", "
    );
  }
}

## Finally run sleuth prep
so <- base::eval(
  base::parse(
    text = command
  )
);


base::saveRDS(
  object = so,
  file = snakemake@output[["rds"]]
);
