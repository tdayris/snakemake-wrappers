#!/usr/bin/R

# This script takes a deseq2 dataset object and performs
# a negative binomial wald test on it

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

# Differential Gene expression
base::library(package = "DESeq2", quietly = TRUE);

# Load DESeq2 dataset
dds_path <- base::as.character(
  x = snakemake@input[["dds"]]
);
dds <- base::readRDS(file = dds_path);

# Build extra parameters for DESeq2 nbinomWaldTest
extra <- "";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["extra"]])
  );
}
base::message("Libraries and dataset loaded");

# Create object
wald <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::nbinomWaldTest(object = dds", extra, ");"
    )
  )
);

# Save results
output_rds <- base::as.character(
  x = snakemake@output[["rds"]]
);

base::saveRDS(
  obj = wald,
  file = output_rds
);
base::message("Wald test over, RDS saved");


names <- DESeq2::resultsNames(
  object = wald
);

output_prefix <- snakemake@output[["tsv"]];
if (! base::file.exists(output_prefix)) {
  base::dir.create(
    path = output_prefix,
    recursive = TRUE
  );
}

alpha_threshold <- 0.5;
if ("alpha_threshold" %in% names(snakemake@params)) {
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

for (resultname in names) {
  base::message(base::paste("Saving results for", resultname))
  results_frame <- DESeq2::results(
    object = wald,
    name = resultname,
    independentFiltering = TRUE,
    # alpha = alpha_threshold,
    # lfcThreshold = fc_threshold,
    pAdjustMethod = "BH",
    cooksCutoff = TRUE
  );

  results_path <- base::file.path(
    output_prefix,
    base::paste0("Deseq2_", resultname, ".tsv")
  );

  utils::write.table(
    x = results_frame,
    file = results_path,
    quote = FALSE,
    sep = "\t",
    row.names = TRUE
  );
}

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
