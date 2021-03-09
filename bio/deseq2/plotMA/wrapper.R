#!/usr/bin/R

# This script takes a deseq2 dataset object and performs
# a mean average plot

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

# Differential Gene expression
base::library(package = "DESeq2", quietly = TRUE);
# Graphic library
base::library(package = "Cairo", quietly = TRUE);

# Load tsv file
res <- utils::read.table(
  file = snakemake@input[["res"]],
  sep = "\t",
  stringsAsFactors = FALSE
);

alpha_threshold <- 0.05;
if ("alpha_threshold" %in% names(snakemake@params)) {
  alpha_threshold <- base::as.numeric(
    x = snakemake@params[["alpha_threshold"]]
  );
}

res$Sig <- res$padj < alpha_threshold;
res <- res[, c("baseMean", "log2FoldChange", "Sig")];

# Build extra parameters for DESeq2 plotMA
extra <- "res";
if ("extra" %in% snakemake@params) {
  extra <- base::paste(
    extra,
    base::as.character(x = snakemake@params[["extra"]]),
    sep = ", "
  );
}

command <- base::paste0(
  "DESeq2::plotMA(",
  extra,
  ");"
);

base::message(command);

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

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
