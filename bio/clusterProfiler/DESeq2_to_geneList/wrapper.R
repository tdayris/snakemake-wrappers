#!/usr/bin/R

# This script takes a deseq2 tsv result and build
# a clusterProfiler compatible geneList object

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");


# Loading DESeq2 for dds object handling
base::library(package = "DESeq2", quietly = TRUE);
# Loading databases
base::library(package = "org.Hs.eg.db", quietly = TRUE);
base::library(package = "org.Mm.eg.db", quietly = TRUE);

# Loading input file
rds <- base::readRDS(
  file=snakemake@input[["rds"]]
);

organism <- org.Hs.eg.db;
if ("organism" %in% base::names(snakemake@params)) {
  if (snakemake@params[["organism"]] == "Mm") {
    organism <- org.Mm.eg.db;
  }
}

# Setting alpha and fold change thresholds
alpha_threshold <- 0.05
if ("alpha_threshold" %in% base::names(snakemake@params)) {
  alpha_threshold <- base::as.numeric(
      x = snakemake@params[["alpha_threshold"]]
  );
}
fc_threshold <- 0.001;
if ("fc_threshold" %in% names(snakemake@params)) {
  fc_threshold <- base::as.numeric(
    x = snakemake@params[["fc_threshold"]]
  );
}

gene_id <- "ENSEMBL";
if ("gene_id" %in% base::names(snakemake@params)) {
  gene_id <- base::as.character(
    x = snakemake@params[["gene_id"]]
  );
}
base::message("Dataset and libraries loaded");

# Gathering results contained within the object
res_names <- DESeq2::resultsNames(
  object = rds
);

# Building output directory
if (! base::file.exists(snakemake@output[["gene_lists"]])) {
  base::dir.create(snakemake@output[["gene_lists"]]);
}

# Building geneLists iteratively
for (resultname in res_names) {
  base::message(base::paste("Building geneList for", resultname))

  out_name <- base::file.path(
    snakemake@output[["gene_lists"]],
    base::paste0(resultname, ".RDS")
  );

  res <- DESeq2::results(
    object = rds,
    name = resultname,
    independentFiltering = TRUE,
    alpha = alpha_threshold,
    lfcThreshold = fc_threshold,
    pAdjustMethod = "BH",
    cooksCutoff = TRUE
  );

  # Adding ENTREZ identifiers
  res$entrez <- mapIds(
    organism,
    keys=row.names(rds),
    column="ENTREZID",
    keytype=gene_id,
    multiVals="first"
  );

  # Buiding geneList object
  geneList <- res[, "log2FoldChange"];
  base::names(geneList) <- res$entrez;
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
    file = out_name
  );
}
base::message("Process over");

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
