#!/bin/R

# Snakemake wrapper for sleuth lrt
# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

base::library(package = "sleuth", quietly = TRUE);

so <- base::readRDS(
    file = snakemake@input[["rds"]]
);

extra <- "";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::as.character(
    x = snakemake@params[["extra"]]
  );
}


# For each beta, run sleuth_wt
for (beta in colnames(so$design_matrix)) {

  # Update extra parameters
  extra_beta <- base::paste0(
    extra,
    ", which_beta = ",
    beta
  );

  # Build command line
  command <- base::paste(
    "sleuth::sleuth_wt(",
    "obj = so",
    extra_beta,
    ")",
    sep = ", "
  );

  # Execute wald test
  so_wt <- base::eval(
    base::parse(
      text = command
    )
  );
}

# Save results
base::saveRDS(
  object = so_wt,
  file = snakemake@output[["rds"]]
);
