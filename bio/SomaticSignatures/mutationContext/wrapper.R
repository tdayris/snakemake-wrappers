#!/usr/bin/R

# Read a Vranged-VCF and a parsed fasta file to build mutation context
library(package = "Rsamtools", quietly = TRUE);
library(package = "VariantAnnotation", quietly = TRUE);
library(package = "SomaticSignatures", quietly = TRUE);

# Load input datasets
fasta <- base::readRDS(
  file = base::as.character(x = snakemake@input[["sequence"]])
);

calling <- base::readRDS(
  file = base::as.character(x = snakemake@input[["call"]])
);

extra <- "vr = calling, ref = fasta";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
}

command <- base::paste0(
  "SomaticSignatures::mutationContext(",
  extra,
  ")"
);
print(command);


# Build context
context <- base::eval(
  base::parse(
    text = command
  )
);

# Filtering empty genotypes
context <- context[context$GT != "./.", ];

# Saving results on user's demand
if ("rds" %in% base::names(snakemake@output)) {
  base::saveRDS(
    object = context,
    file = base::as.character(x = snakemake@output[["rds"]])
  );
}

if ("tsv" %in% base::names(snakemake@output)) {
  utils::write.table(
    x = base::as.data.frame(context),
    file = base::as.character(x = snakemake@output[["tsv"]]),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  );
}
