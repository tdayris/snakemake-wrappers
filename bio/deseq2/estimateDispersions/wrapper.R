#!/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Perform sample dispertion estimates

# Differential Gene expression
base::library(package = "DESeq2", quietly = TRUE);

# Cast input path as character
dds_path <- base::as.character(x = snakemake@input[["dds"]]);
dds <- base::readRDS(dds_path);

# Cast locfunc as function name
extra <- "";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["extra"]])
  );
}
# Create object
dds <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::estimateDispersions(dds", extra, ");"
    )
  )
);

# Save as RDS
output_path <- base::as.character(x = snakemake@output[["disp"]]);
base::saveRDS(
  obj = dds,
  file = output_path
);
