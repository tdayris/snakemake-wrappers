#!/usr/bin/R

# This script takes a deseq2 dataset object and performs
# a default DESeq2 analysis

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Differential Gene expression
base::library(package = "SummarizedExperiment", quietly = TRUE);
base::library(package = "DESeq2", quietly = TRUE);

cleanColData <- function(dds, factor) {
  NApos <- base::is.na(dds[[factor]]);
  dds2 <- dds[, !NApos];
  for (coldata in base::names(colData(dds2))) {
    dds2[[coldata]] <- dds[[coldata]][!NApos]
  }
  return(dds2)
}

# Load DESeq2 dataset
dds_path <- base::as.character(x = snakemake@input[["dds"]]);
dds <- base::readRDS(file = dds_path);

# Build extra parameters for DESeq2
extra_deseq2 <- "object = dds";
if ("extra" %in% snakemake@params) {
  extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["extra"]])
  );
}
base::message("Libraries and dataset loaded");

deseq2_cmd <- base::paste0(
  "DESeq2::DESeq(", extra_deseq2, ");"
);
base::message("DESeq2 command line:");
base::message(deseq2_cmd);

# Create object
wald <- tryCatch({
    base::eval(base::parse(text = deseq2_cmd))
  },
  error = function(e) {
    dds <- cleanColData(
      dds = dds,
      factor = base::as.character(x = snakemake@params[["factor"]])
    );
    base::eval(base::parse(text = deseq2_cmd))
  }
);
#wald <- base::eval(base::parse(text = deseq2_cmd));

# Save results as RDS
output_rds <- base::as.character(x = snakemake@output[["rds"]]);
base::saveRDS(obj = wald, file = output_rds);
base::message("Wald test over, RDS saved");

# Saving results as TSV
names <- DESeq2::resultsNames(object = wald);

if ("deseq2_result_dir" %in% base::names(snakemake@output)) {
  # Recovreing extra parameters for TSV tables
  extra_results <- "object = wald, name = resultname";
  if ("extra_results" %in% base::names(snakemake@params)) {
    extra_results <- base::paste(
      ", ",
      base::as.character(x = snakemake@params[["extra_results"]])
    );
  }

  # DESeq2 result dir will contain all results available in the Wald object
  output_prefix <- snakemake@output[["deseq2_result_dir"]];
  if (! base::file.exists(output_prefix)) {
    base::dir.create(path = output_prefix,recursive = TRUE);
  }

  results_cmd <- base::paste0("DESeq2::results(", extra_results, ")");
  base::message("Command line used for TSV results creation:");
  base::message(results_cmd);

  for (resultname in names) {
    # Building table
    base::message(base::paste("Saving results for", resultname))
    results_frame <- base::eval(base::parse(text = results_cmd));

    results_path <- base::file.path(
      output_prefix,
      base::paste0("Deseq2_", resultname, ".tsv")
    );

    # Saving table
    utils::write.table(
      x = results_frame,
      file = results_path,
      quote = FALSE,
      sep = "\t",
      row.names = TRUE
    );
  }
}

if ("deseq2_tsv" %in% base::names(snakemake@output)) {
  base::message("Saving results as TSV");
  # This will save results based on their contrast identifier
  # For more information, see deseq2 vignette or type ?results
  constrast <- c(
    base::as.character(x = snakemake@params[["factor"]]),
    base::as.character(x = snakemake@params[["numerator"]]),
    base::as.character(x = snakemake@params[["denominator"]])
  );
  extra_results <- base::paste(
    "object=wald",
    "contrast=constrast",
    sep=", "
  );
  if ("extra_results" %in% base::names(snakemake@params)) {
    extra_results <- base::paste(
      extra_results,
      base::as.character(x = snakemake@params[["extra_results"]]),
      sep=", "
    );
  }
  results_cmd <- base::paste0("DESeq2::results(", extra_results, ")");
  base::message("Result extraction command:");
  results_frame <- base::eval(base::parse(text = results_cmd));

  # Saving table
  utils::write.table(
    x = results_frame,
    file = base::as.character(x = snakemake@output[["deseq2_tsv"]]),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE
  );
}

# Saving normalized counts on demand
table <- SummarizedExperiment::assay(wald);
if ("normalized_counts" %in% base::names(snakemake@output)) {
  output_table <- base::as.character(x=snakemake@output[["normalized_counts"]]);
  utils::write.table(x = table, file = output_table, sep = "\t", quote = FALSE);
  base::message("Normalized counts saved as TSV");
}

if ("dst" %in% base::names(snakemake@output)) {
  output_rds <- base::as.character(
    x=snakemake@output[["dst"]]
  );
  base::saveRDS(
    obj = table,
    file = output_rds
  );
}
