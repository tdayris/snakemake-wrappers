#/usr/bin/env R

# Snakemake wrapper for biological translator

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open = "wt")
base::sink(log.file)
base::sink(log.file, type = "message")

base::library(package = "AnnotationDbi", character.only = TRUE)
base::library(package = "org.Hs.eg.db", character.only = TRUE)
base::library(package = "org.Mm.eg.db", character.only = TRUE)

orgdb <- "org.Hs.eg.db"
if ("orgdb" %in% base::names(snakemake@params)) {
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

save_tsv <- function(data, keytype, weights, out_key) {
    out_path <- base::as.character(x = snakemake@out_key)
    base::message("Saving TSV to ", out_path)

    utils::write.table(
        x = data[c(keytype, weights), ]
        file = out_path,
        sep = "\t"
    )
}

save_rds <- function(data, keytype, weights, out_key) {
  out_path <- base::as.character(x = snakemake@out_key)
  base::message("Saving RDS to ", out_path)

  geneList <- build_gene_list(
    gene_frame = data,
    comparison_name = weights,
    keyid = keytype
  )

  base::saveRDS(
    object = geneList,
    file = out_path
  )
}


save_universe <- function(data, keytype, out_key) {
  out_path <- base::as.character(x = snakemake@out_key)
  base::message("Saving RDS to ", out_path)

  base::saveRDS(
    object = data[keytype, ],
    file = out_path
  )
}


# Loading input file
tsv <- utils::read.table(
  file = base::as.character(x = snakemake@input[["tsv"]]),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
base::message("Dataset and libraries loaded")

gene_id_type <- get_parameter("gene_id_type", "ENSEMBL")
base::message("KeyType acquired")

if (! "ENSEMBL" %in% colnames(tsv)) {
  base::message("Adding ENSEMBL keys to the gene table")
  tsv$ENSEMBL <- mapIds(
    orgdb,
    keys = tsv$SYMBOL,
    colum = 'ENSEMBL',
    keytype = "SYMBOL",
    multiVals = 'first'
  )
}

if (! "SYMBOL" %in% colnames(tsv)) {
  base::message("Adding SYMBOL keys to the gene table")
  tsv$SYMBOL <- mapIds(
    orgdb,
    keys = tsv$ENSEMBL,
    colum = 'SYMBOL',
    keytype = "ENSEMBL",
    multiVals = 'first'
  )
}

base::message("Adding ENTREZID keys to the gene table ", gene_id_type)
tsv$ENTREZID <- mapIds(
  orgdb,
  keys = tsv$ENSEMBL,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

base::message("Adding ENSEMBLPROT keys to the gene table")
tsv$ENSEMBLPROT <- mapIds(
  orgdb,
  keys = tsv$ENSEMBL,
  column = "ENSEMBLPROT",
  keytype = "ENSEMBL",
  multiVals = "first"
)
base::print(utils::head(tsv))

weights <- base::colnames(tsv)
weights <- comparison_names[! comparison_names %in% c("SYMBOL", "ENSEMBL", "ENTREZID", "ENSEMBLPROT")]

if ("tsv_ensembl" %in% base::names(snakemake@output)) {
    save_tsv <- function(tsv, "ENSEMBL", weights, "tsv_ensembl")
}

if ("tsv_entrez" %in% base::names(snakemake@output)) {
    save_tsv <- function(tsv, "ENTREZID", weights, "tsv_entrez")
}

if ("tsv_ensemblprot" %in% base::names(snakemake@output)) {
    save_tsv <- function(tsv, "ENSEMBLPROT", weights, "tsv_ensemblprot")
}

if ("tsv_symbol" %in% base::names(snakemake@output)) {
    save_tsv <- function(tsv, "SYMBOL", weights, "tsv_symbol")
}

if ("rds_ensembl" %in% base::names(snakemake@output)) {
    save_rds <- function(tsv, "ENSEMBL", weights, "rds_ensembl")
}

if ("rds_entrez" %in% base::names(snakemake@output)) {
    save_rds <- function(tsv, "ENTREZID", weights, "rds_entrez")
}

if ("rds_ensemblprot" %in% base::names(snakemake@output)) {
    save_rds <- function(tsv, "ENSEMBLPROT", weights, "rds_ensemblprot")
}

if ("rds_symbol" %in% base::names(snakemake@output)) {
    save_rds <- function(tsv, "SYMBOL", weights, "rds_symbol")
}

if ("universe_symbol" %in% base::names(snakemake@output)) {
    save_universe <- function(tsv, "SYMBOL", "universe_symbol")
}

if ("universe_ensembl" %in% base::names(snakemake@output)) {
    save_universe <- function(tsv, "ENSEMBL", "universe_ensembl")
}

if ("universe_ensemblprot" %in% base::names(snakemake@output)) {
    save_universe <- function(tsv, "ENSEMBLPROT", "universe_ensemblprot")
}

if ("universe_entrez" %in% base::names(snakemake@output)) {
    save_universe <- function(tsv, "ENTREZID", "universe_entrez")
}


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()
