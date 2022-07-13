#!/usr/bin/R

# This script takes a geneList object and performs
# an enrichment analysis based on msigdb

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

# MSigDB package used to load the pathways
base::library(package = "msigdbr", quietly = TRUE);
# Perform gene enrichment
base::library(package = "clusterProfiler", quietly = TRUE);
# Handling large datasets
base::library(package = "dplyr", quietly = TRUE);
# Loading databases
base::library(package = "org.Hs.eg.db", quietly = TRUE);
base::library(package = "org.Mm.eg.db", quietly = TRUE);


# Loading input dataset
geneList <- base::readRDS(
  file = snakemake@input[["rds"]]
);

organism <- org.Hs.eg.db;
spieces <- "Homo sapiens";
if ("organism" %in% base::names(snakemake@params)) {
  if (snakemake@params[["organism"]] == "Mm") {
    organism <- org.Mm.eg.db;
    spieces = "Mus musculus";
  }
}

# Loading MSigDB
msigdb_extra <- 'species = spieces';
if ("msigdb_extra" %in% base::names(snakemake@params)) {
  msigdb_extra <- base::paste(
    msigdb_extra,
    base::as.character(x = snakemake@params[["msigdb_extra"]]),
    sep=", "
  )
}

command <- base::paste0(
  "msigdbr::msigdbr(",
  msigdb_extra,
  ")"
)
base::message(command);

db <- base::eval(
  base::parse(
    text = command
  )
) %>% dplyr::select(gs_name, entrez_gene)
base::message("Libraries, datasets and database retrieved");

# Performing GSEA
gsea_extra <- "geneList, TERM2GENE = db";
if ("msigdb_gsea_extra" %in% base::names(snakemake@params)) {
  gsea_extra <- base::paste(
    gsea_extra,
    base::as.character(snakemake@params[["msigdb_gsea_extra"]]),
    sep=", "
  );
}
gsea <- base::eval(
  base::parse(
    text = command
  )
)

gsea <- clusterProfiler::setReadable(
    gsea,
    OrgDb = organism,
    keyType = "ENTREZID"
);

# Saving results
base::saveRDS(
  obj = gsea,
  file = snakemake@output[["rds"]]
);

utils::write.table(
  x = gsea,
  file = snakemake@output[["tsv"]]
);

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
