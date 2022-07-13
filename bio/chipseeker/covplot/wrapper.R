#!/usr/bin/R

# This script loads peaks into R with ChipSeeker

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"


# Load libraries
base::library(package = "ChIPseeker", quietly = TRUE);

# Read peaks
tagMatrixList <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
);
base::message("Dataset loaded");

extra <- "peak = tagMatrixList";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
}

# Create R command line
command <- base::eval(
  base::parse(
    text = base::paste0(
      "ChIPseeker::covplot(", extra, ")"
    )
  )
);


# Build plot
png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px",
  type = "cairo"
);

base::eval(
  base::parse(
    text = command
  )
);

dev.off();
