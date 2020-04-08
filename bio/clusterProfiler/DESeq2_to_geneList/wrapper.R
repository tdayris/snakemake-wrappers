#!/usr/bin/R

# This script takes a deseq2 tsv result and build
# a clusterProfiler compatible geneList object

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"


# Loading databases
base::library(package = "org.Hs.eg.db", quietly = TRUE);
base::library(package = "org.Mm.eg.db", quietly = TRUE);

# Loading input file
res <- utils::read.table(
  file = snakemake@input[["res"]],
  sep = "\t",
  header = TRUE
);


# Setting alpha thresholds
alpha_threshold <- 0.05
if ("alpha_threshold" %in% base::names(snakemake@params)) {
  alpha_threshold <- base::as.numeric(
      x = snakemake@params[["alpha_threshold"]]
  )
}

# Subsetting initial data
res <- res[res$padj <= alpha_threshold, ];

# Buiding geneList object
geneList <- res[, "log2FoldChange"];
base::names(geneList) <- base::names(res);
geneList <- sort(
  x = geneList,
  decreasing = TRUE,
  na.last = NA
);

base::message(
  utils::head(x = geneList)
);

base::saveRDS(
  object = geneList,
  file = snakemake@output[["rds"]]
);
