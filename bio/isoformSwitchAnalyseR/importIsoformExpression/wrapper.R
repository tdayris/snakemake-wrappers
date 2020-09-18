#!/usr/bin/R

# This script takes an list of Salmon directories
# and builds a IsoformSwitchAnalyseR object

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"


# The main package
base::library(package = "IsoformSwitchAnalyzeR", quietly = TRUE);

# Gathering input files
quant_files <- sapply(
  snakemake@input[["quant"]],
  function(quant) base::as.character(x = quant)
);

design_path <- base::as.character(
  x = snakemake@input[["design"]]
);

separator <- base::ifelse(
  base::endsWith(design_path, ".csv"),
  ",",
  "\t"
);

design_table <- utils::read.table(
  file=design_path,
  sep=separator,
  header=TRUE
);

sample_id_col <- "Sample_id";
if ("sample_id_col" %in% names(snakemake@params)) {
  sample_id_col <- snakemake@params[["sample_id_col"]];
}

condition_col <- "condition";
if ("condition" %in% names(snakemake@params)) {
  condition_col <- snakemake@params[["condition_col"]];
}

design <- data.frame(
  sampleID = design_table[, sample_id_col],
  condition = design_table[, condition_col]
);

gtf_path <- base::as.character(x = snakemake@input[["gtf"]]);
fasta_path <- base::as.character(x = snakemake@input[["fasta"]]);

# Guessing the argument type in importIsoformExpression
import_type <- ifelse(
  all(file_test("-f", quant_files)),
  "sampleVector",
  "parentDir"
);


# Building command line arguments
extra <- base::paste0(import_type, " = quant_files");
if ("extra_isoform_expression" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra_isoform_expr"]],
    sep = ", "
  );
}

# Buiding command line itself
command <- base::paste0(
  "IsoformSwitchAnalyzeR::importIsoformExpression(",
  extra,
  ")"
);
base::message("Libraries and input data loaded");
base::message(command);

# Running command
quantification <- base::eval(
  base::parse(
    text = command
  )
);

new_names <- sapply(
  colnames(quantification$counts),
  function(sampleid) base::ifelse(
      sampleid == "isoform_id",
      "isoform_id",
      base::basename(base::dirname(sampleid))
    )
);
colnames(quantification$counts) <- new_names;
colnames(quantification$abundance) <- new_names;

command = base::paste0(
  "IsoformSwitchAnalyzeR::importRdata(",
  "isoformCountMatrix = quantification$counts, ",
  "isoformRepExpression = quantification$abundance, ",
  "designMatrix = design, ",
  "isoformExonAnnoation = gtf_path, ",
  "isoformNtFasta = fasta_path "
);
command = "IsoformSwitchAnalyzeR::importRdata("
if ("extra_rdata" %in% names(snakemake@params)) {
  command = base::paste0(command, ", ", snakemake@params[["extra_rdata"]], ")");
} else {
  command = base::paste0(command, ")");
}

base::message("Import R data:");
base::message(command);

# Running command
switch_list <- base::eval(
  base::parse(
    text = command
  )
);

# switch_list <- IsoformSwitchAnalyzeR::importRdata(
#   isoformCountMatrix = quantification$counts,
#   isoformRepExpression = quantification$abundance,
#   designMatrix = design,
#   isoformExonAnnoation = gtf_path,
#   isoformNtFasta = fasta_path,
#   showProgress = TRUE
# );

print("Switchlist built")

# Saving results
base::saveRDS(
  switch_list,
  file = snakemake@output[["rds"]]
);
