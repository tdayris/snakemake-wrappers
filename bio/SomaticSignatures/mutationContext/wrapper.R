#!/usr/bin/R

# Read a Vranged-VCF and a parsed fasta file to build mutation context
library(package = "Rsamtools", quietly = TRUE);
library(package = "VariantAnnotation", quietly = TRUE);
library(package = "SomaticSignatures", quietly = TRUE);
print("loaded!")

# Load input datasets
fasta <- base::readRDS(
  file = base::as.character(x = snakemake@input[["sequence"]])
);
print(fasta)

calling <- base::readRDS(
  file = base::as.character(x = snakemake@input[["call"]])
);
print(calling)

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

base::saveRDS(
  object = context,
  file = base::as.character(x = snakemake@output[["rds"]])
);
