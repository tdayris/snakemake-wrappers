#!/usr/bin/R

# Load VCF file and save a Range object RDS-formatted.
base::library(package="VariantAnnotation", quietly=TRUE);

vcf_file <- VariantAnnotation::VcfFile(
  file = base::as.character(x = snakemake@input[["calls"]]),
  index = base::as.character(x = snakemake@input[["tbi"]])
);

print(vcf_file)

extra <- "x = vcf_file";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
    extra,
    base::as.character(x = snakemake@params["extra"]),
    sep=", "
  );
}

command <- base::paste0(
  "VariantAnnotation::readVcfAsVRanges(",
  extra,
  ");"
);
print(command)

# Load vcf file as Range object
vcf <- base::eval(
  base::parse(
    text = command
  )
);


# Save as RDS
base::saveRDS(
  object = vcf,
  file = base::as.character(snakemake@output[["rds"]])
);
