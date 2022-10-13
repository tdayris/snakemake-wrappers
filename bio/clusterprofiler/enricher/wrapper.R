#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2022, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for ClusterProfiler enrichment analysis

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]], open = "wt")
base::sink(log.file)
base::sink(log.file, type = "message")

# Read snakemake input, and guess wether it's a RDS or a TSV
read_input <- function(data) {
    data <- base::as.character(x = data)
    if (base::endsWith(x = data, suffix = ".RDS") {
        data <- readRDS(file = data)
    } else {
        data <- utils::read.table(
            file = data,
            header = FALSE,
            sep = "\t",
            stringsAsFactors = FALSE
        )
    }
    base::return(data)
}

# Load gene id list and optional fold changes
gene <- NULL
fold_change <- NULL
if ("fold_change" %in% base::names(x = snakemake@input)) {
    # Then user provided a fold_change file. 
    # Heatmap will be coloured, and input.gene
    # becomes optional.
    fold_change <- read_input(data = snakemake@input[["fold_change"]])

    if ("gene" %in% base::names(x = snakemake@input)) {
        gene <- read_input(data = snakemake@input[["gene"]])
    } else {
        gene <- fold_change[0, ]
    }
} else {
    # Then user did not provide any fold_change file.
    # input.gene is required, and no colors will be
    # available on heatmaps.
}