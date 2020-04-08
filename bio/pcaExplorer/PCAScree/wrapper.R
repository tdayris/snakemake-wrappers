#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script takes a deseq2 transform object and performs
# a plot over pca loadings


base::library(package = "DESeq2");        # Differential analysis
base::library(package = "pcaExplorer");   # Handles PCAs
base::library(package = "Cairo");         # Graphic library

# Overload output defaults in order to avoid
# X11 foreward errors on cluster nodes
options(bitmapType="cairo");

# Load specified input files
dst_path <- base::as.character(
  x = snakemake@input[["dst"]]
);
dst <- base::readRDS(file = dst_path);
pca <- stats::prcomp(t(SummarizedExperiment::assay(dst)))


# Load extra parameters
extra <- "pca"
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
      extra,
      snakemake@params[["extra"]],
      sep = ", "
  );
}

# Build plot
png(
  filename = snakemake@output[["png"]],
  width = 1024,
  height = 768,
  units = "px",
  type = "cairo"
);

command <- base::paste0(
  "pcaExplorer::pcascree(",
  extra,
  ");"
);

base::message(command);

base::eval(
  base::parse(
    text = command
  )
);

dev.off()
