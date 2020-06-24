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


# For each beta, run sleuth_wt
for (beta in colnames(so$design_matrix)) {

  # Update extra parameters
  extra_beta <- base::paste0(
    extra,
    ", test = ",
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

  # Extract results
  so_results <- base::eval(
    base::parse(
      text = command
    )
  );

  # Save table
  if ("gene_analysis" %in% names(snakemake@params)) {
    
  }
}
