#/usr/bin/env R

# Snakemake wrapper for biological translator

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

base::library(package = "AnnotationDbi", quietly = TRUE);
base::library(package="org.Hs.eg.db", quietly=TRUE);
base::library(package="org.Mm.eg.db", quietly=TRUE);
base::library(package = "clusterProfiler", quietly = TRUE);

orgdb <- "org.Hs.eg.db";
if ( "orgdb" %in% base::names(snakemake@params)) {
  orgdb <- base::as.character(x = snakemake@params[["orgdb"]])
}

get_parameter <- function(param_name, default_value) {
  # Return the provided parameter or a default value otherwise
  res <- default_value;
  if (param_name %in% base::names(snakemake@params)) {
    res <- base::as.character(x = snakemake@params[[param_name]]);
  }
  return(res)
}

build_gene_list <- function(gene_frame, comparison_name, keyid) {
  tmp <- gene_frame[, c(comparison_name, keyid)];
  tmp <- tmp[! is.na(tmp[, comparison_name]), ];
  tmp <- tmp[! is.na(tmp[, keyid]), ];

  geneList <- tmp[, comparison_name];
  base::names(geneList) <- tmp[, keyid];
  geneList <- geneList[order(geneList, decreasing=TRUE)];
  return(geneList);
}


# Loading input file
tsv <- utils::read.table(
  file = base::as.character(x = snakemake@input[["tsv"]]),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
);

gene_id_type <- get_parameter("gene_id_type", "ENSEMBL");
base::message("Dataset and libraries loaded");

if (! "ENSEMBL" %in% colnames(tsv)) {
  base::message("Adding ENSEMBL keys to the gene table");
  tsv$ENSEMBL <- mapIds(
    org.Hs.eg.db,
    keys=tsv$SYMBOL,
    colum='ENSEMBL',
    keytype="SYMBOL",
    multiVals='first'
  );
}

 if (! "SYMBOL" %in% colnames(tsv)) {
  base::message("Adding SYMBOL keys to the gene table");
  tsv$SYMBOL <- mapIds(
    org.Hs.eg.db,
    keys=tsv$ENSEMBL,
    colum='SYMBOL',
    keytype="ENSEMBL",
    multiVals='first'
  );
}

base::message("Adding ENTREZID keys to the gene table ", gene_id_type);
tsv$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys=tsv$ENSEMBL,
  column="ENTREZID",
  keytype="ENSEMBL",
  multiVals="first"
);

base::message("Adding ENSEMBLPROT keys to the gene table");
tsv$ENSEMBLPROT <- mapIds(
  org.Hs.eg.db,
  keys=tsv$ENSEMBL,
  column="ENSEMBLPROT",
  keytype="ENSEMBL",
  multiVals="first"
);
print(head(tsv))

comparison_names <- base::colnames(tsv);
comparison_names <- comparison_names[! comparison_names %in% c("SYMBOL", "ENSEMBL", "ENTREZID", "ENSEMBLPROT")];
nb_comparisons <- base::length(comparison_names);
base::message("Comparison acquired, ", comparison_names);

for (factor in 1:nb_comparisons) {
  base::message("Working on comparison:", comparison_names[factor])

  if ("tsv" %in% base::names(snakemake@output)) {
    file_path <- base::as.character(
      snakemake@output[["tsv"]][factor]
    );
    base::message("Saving TSV to ", file_path);

    geneList <- build_gene_list(
      gene_frame=tsv,
      comparison_name=comparison_names[factor],
      keyid="SYMBOL"
    );

    utils::write.table(
      x=as.data.frame(geneList),
      file=file_path,
      sep="\t"
    );
  }

  if ("symbol_rds" %in% base::names(snakemake@output)) {
    rds_path <- base::as.character(
      snakemake@output[["symbol_rds"]][factor]
    );
    base::message("Saving SYMBOLS to ", rds_path);

    geneList <- build_gene_list(
      gene_frame=tsv,
      comparison_name=comparison_names[factor],
      keyid="SYMBOL"
    );

    base::saveRDS(
      object=geneList,
      file=rds_path
    );
  }

  if ("ensembl_rds" %in% base::names(snakemake@output)) {
    rds_path <- base::as.character(
      snakemake@output[["ensembl_rds"]][factor]
    );
    base::message("Saving ENSEMBL to ", rds_path);

    geneList <- build_gene_list(
      gene_frame=tsv,
      comparison_name=comparison_names[factor],
      keyid="ENSEMBL"
    );

    base::saveRDS(
      object=geneList,
      file=rds_path
    );
  }

  if ("protein_rds" %in% base::names(snakemake@output)) {
    rds_path <- base::as.character(
      x=snakemake@output[["protein_rds"]][factor]
    );
    base::message("Saving Protein RDS to ", rds_path);

    proteinList <- build_gene_list(
      gene_frame=tsv,
      comparison_name=comparison_names[factor],
      keyid="ENSEMBLPROT"
    );

    base::saveRDS(
      object=proteinList,
      file=rds_path
    );
  }

  if ("entrez_rds" %in% base::names(snakemake@output)) {
    rds_path <- base::as.character(
      x=snakemake@output[["entrez_rds"]][factor]
    );
    base::message("Saving Gene RDS to ", rds_path);

    geneList <- build_gene_list(
      gene_frame=tsv,
      comparison_name=comparison_names[factor],
      keyid="ENTREZID"
    );

    base::saveRDS(
      object=geneList,
      file=rds_path
    );
  }

  if ("universe" %in% base::names(snakemake@output)) {
    rds_path <- base::as.character(
      x=snakemake@output[["universe"]][factor]
    );
    base::message("Saving universe to ", rds_path);

    geneList <- tsv$ENTREZID;
    names(geneList) <- tsv$SYMBOL;

    base::saveRDS(
      object=geneList,
      file=rds_path
    );
  }
}

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
