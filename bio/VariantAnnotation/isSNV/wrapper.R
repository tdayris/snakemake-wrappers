#!/usr/bin/R

# Filter VCF based on the SNP status
library(package = "VariantAnnotation", quietly = TRUE);
library(package = "Biostrings", quietly = TRUE);

# Load dataset
vcf <- base::readRDS(
  base::as.character(x = snakemake@input[["calls"]])
);

# Filter dataset on SNVs
keep_snv <- VariantAnnotation::isSNV(vcf, singleAltOnly=TRUE);
keep_ref <- VariantAnnotation::ref(vcf) %in% Biostrings::DNA_BASES;
keep_alt <- VariantAnnotation::alt(vcf) %in% Biostrings::DNA_BASES;
vcf <- vcf[keep_alt & keep_ref & keep_snv, ];

# Save results
base::saveRDS(
  object = vcf,
  file = base::as.character(x = snakemake@output[["rds"]])
);
