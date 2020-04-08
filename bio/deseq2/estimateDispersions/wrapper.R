#!/bin/R
# among samples for further normalization

base::library("DESeq2");     # Differential Gene expression

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
