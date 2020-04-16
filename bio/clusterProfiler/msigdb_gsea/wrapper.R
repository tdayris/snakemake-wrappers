#!/usr/bin/R

# This script takes a geneList object and performs
# an enrichment analysis based on msigdb

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

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
if ("organism" %in% base::names(snakemake@params)) {
  if (snakemake@params[["organism"]] == "Mm") {
    organism <- org.Mm.eg.db;
  }
}

# Loading MSigDB
msigdb_extra <- 'species = "Homo sapiens"';
if ("msigdb_extra" %in% snakemake@params) {
  msigdb_extra <- base::as.character(
    x = snakemake@params[["msigdb_extra"]]
  );
}

command <- base::paste0(
  "msigdbr::msigdbr(",
  extra,
  ")"
)

db <- base::eval(
  base::parse(
    text = command
  )
) %>% dplyr::select(gs_name, entrez_gene)
base::message("Libraries, datasets and database retrieved");

# Performing GSEA
gsea <- clusterProfiler::GSEA(
  geneList,
  TERM2GENE = db
);

gsea <- clusterProfiler::setReadable(
    gsea,
    organism,
    keytype = "ENTREZ"
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
