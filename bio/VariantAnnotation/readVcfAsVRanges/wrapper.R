#!/usr/bin/R

# Load VCF file and save a Range object RDS-formatted.
base::library(package="VariantAnnotation", quietly=TRUE);

calling <- VariantAnnotation::VcfFile(
  file = base::as.character(x = snakemake@input[["calls"]]),
  index = base::as.character(x = snakemake@input[["tbi"]])
);

extra <- "x = calling";
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
vrange <- base::eval(
  base::parse(
    text = command
  )
);


# Save as RDS
base::saveRDS(
  object = vrange,
  file = base::as.character(snakemake@output[["rds"]])
);
