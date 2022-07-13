#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script takes a deseq2 dataset object and performs
# a rlog transformation on it

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

# Differential Gene expression
base::library(package = "DESeq2", quietly = TRUE);

# Cast input path as character
dds_path <- base::as.character(x = snakemake@input[["dds"]]);
dds <- base::readRDS(file = dds_path);

# Recover extra parameters
extra <- "";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
    "object = dds",
    base::as.character(x = snakemake@params[["extra"]]),
    sep = ", "
  );
}

# Create object
rlog <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::rlogTransformation(", extra, ");"
    )
  )
);

tsv <- SummarizedExperiment::assay(rlog);

# Save results
output_rds <- base::as.character(snakemake@output[["rds"]]);
base::saveRDS(
  obj = rlog,
  file = output_rds
);


output_table <- base::as.character(snakemake@output[["tsv"]]);
utils::write.table(
  x = tsv,
  file = output_table,
  sep = "\t",
  quote = FALSE
);

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
