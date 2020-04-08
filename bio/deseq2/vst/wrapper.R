# This script takes a deseq2 dataset object and performs
# a variance stabilizing transformations transformation on it

base::library("DESeq2");                 # Differential Gene expression
base::library("SummarizedExperiment");   # Handle large datasets

# Cast input path as character
dds_path <- base::as.character(x = snakemake@input[["dds"]]);
dds <- base::readRDS(file = dds_path);

# Recover extra parameters
extra <- "";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["extra"]])
  );
}

# Create object
vst <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::vst(object = dds ", extra, ");"
    )
  )
);

# Save results
output_rds <- base::as.character(snakemake@output[["rds"]]);
base::saveRDS(
  obj = vst,
  file = output_rds
);


output_table <- base::as.character(snakemake@output[["tsv"]]);
tsv <- SummarizedExperiment::assay(vst);
utils::write.table(
  x = tsv,
  file = output_table,
  sep = "\t",
  quote = FALSE
);
