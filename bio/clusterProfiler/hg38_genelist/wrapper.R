#/usr/bin/env R

# Snakemake wrapper for biological translator

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

base::library(package = "AnnotationDbi", quietly = TRUE);
base::library(package="org.Hs.eg.db", quietly=TRUE);
base::library(package = "clusterProfiler", quietly = TRUE);
orgdb <- "org.Hs.eg.db";

get_parameter <- function(param_name, default_value) {
  # Return the provided parameter or a default value otherwise
  res <- default_value;
  if (param_name %in% base::names(snakemake@params)) {
    res <- base::as.character(x = snakemake@params[[param_name]]);
  }
  return(res)
}

build_gene_list <- function(gene_frame, comparison_name) {
  geneList <- gene_frame[, comparison_name];
  base::names(geneList) <- gene_frame$ENTREZID;
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

tsv$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys=tsv[, gene_id_type],
  column="ENTREZID",
  keytype=gene_id_type,
  multiVals="first"
)

comparison_names <- base::colnames(tsv);
comparison_names <- comparison_names[! comparison_names %in% c("ENSEMBL", "ENTREZID")];
nb_comparisons <- base::length(comparison_names);
base::message("Comparison acquired");

for (factor in 1:nb_comparisons) {
  base::message("Working on comparison:", comparison_names[factor])
  geneList <- build_gene_list(
    gene_frame=tsv,
    comparison_name=comparison_names[factor]
  );

  if ("tsv" %in% base::names(snakemake@output)) {
    file_path <- base::as.character(
      snakemake@output[["tsv"]][factor]
    );
    base::message("Saving TSV to ", file_path);

    utils::write.table(
      x=as.data.frame(geneList),
      file=file_path,
      sep="\t"
    );
  }

  if ("rds" %in% base::names(snakemake@output)) {
    rds_path <- base::as.character(
      snakemake@output[["rds"]][factor]
    );
    base::message("Saving RDS to ", rds_path);
    print(head(geneList))
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
