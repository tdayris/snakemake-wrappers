#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script takes a deseq2 transform object and performs
# a plot over sample pairwise correlation


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

gene_number <- 100;
if ("gene_number" %in% names(snakemake@params)) {
  gene_number <- base::as.numeric(
    snakemake@params[["gene_number"]]
  );
}


# Load extra parameters
extra <- "counts(dst)[1:gene_number, ]"
if ("extra" %in% names(snakemake@params)) {
  if (snakemake@params[["extra"]] != "") {
    extra <- base::paste(
        extra,
        snakemake@params[["extra"]],
        sep = ", "
    );
  }
}
message(gene_number)

# Build plot
w <- 1024;
if ("w" in base::names(snakemake@params)) {
  w <- base::as.numeric(snakemake@params[["w"]]);
}
h <- 768;
if ("h" in base::names(snakemake@params)) {
  h <- base::as.numeric(snakemake@params[["h"]]);
}

png(
  filename = snakemake@output[["png"]],
  width = w,
  height = h,
  units = "px",
  type = "cairo"
);

command <- base::paste0(
  "pcaExplorer::pair_corr(",
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
