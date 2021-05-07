#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script takes a deseq2 transform object and performs
# a plots over factor correlation with pca axes

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

base::library(package = "DESeq2");        # Differential analysis
base::library(package = "pcaExplorer");   # Handles PCAs
base::library(package = "Cairo");         # Graphic library

cleanColData <- function(dds, factor) {
  NApos <- base::is.na(dds[[factor]]);
  dds2 <- dds[, !NApos];
  for (coldata in base::names(colData(dds2))) {
    dds2[[coldata]] <- dds[[coldata]][!NApos]
    if (length(levels(dds2[[coldata]])) == 1) {dds2[[coldata]] <- NULL}
  }
  return(dds2)
}

# Overload output defaults in order to avoid
# X11 foreward errors on cluster nodes
options(bitmapType="cairo");

# Load specified input files
dst_path <- base::as.character(
  x = snakemake@input[["dst"]]
);
dst <- base::readRDS(file = dst_path);
pca <- stats::prcomp(t(SummarizedExperiment::assay(dst)))

dds_path <- base::as.character(
  x = snakemake@input[["dds"]]
);
dds <- base::readRDS(file = dds_path);

base::message(colData(dds));
base::message(pca);
corrs_pca <- tryCatch({
    pcaExplorer::correlatePCs(pca, colData(dds))
  },
  error = function(e) {
    dds <- cleanColData(
      dds = dds,
      factor = base::as.character(x = snakemake@params[["factor"]])
    );
    pcaExplorer::correlatePCs(pca, colData(dds))
  }
);

# Load extra parameters
extra <- "corrs_pca"
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
      extra,
      snakemake@params[["extra"]],
      sep = ", "
  );
}

command <- base::paste0(
  "pcaExplorer::plotPCcorrs(",
  extra,
  ");"
);

base::message(command);

# Build plot
w <- 1024;
if ("w" %in% base::names(snakemake@params)) {
  w <- base::as.numeric(snakemake@params[["w"]]);
}
h <- 768;
if ("h" %in% base::names(snakemake@params)) {
  h <- base::as.numeric(snakemake@params[["h"]]);
}

png(
  filename = snakemake@output[["png"]],
  width = w,
  height = h,
  units = "px",
  type = "cairo"
);

tryCatch({
    base::eval(base::parse(text = command))
  },
  error = function(e) {
    dds <- cleanColData(
      dds = dds,
      factor = base::as.character(x = snakemake@params[["factor"]])
    );
    base::eval(base::parse(text = command))
  }
);

dev.off();


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
