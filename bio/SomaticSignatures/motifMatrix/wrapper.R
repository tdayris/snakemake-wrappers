#!/usr/bin/R

# Build mutation motif matrix from a VRange-formatted calling
library(package = "SomaticSignatures", quietly=TRUE);

# Loading input dataset
mutations <- base::readRDS(
  file = base::as.character(
    x = snakemake@input[['context']]
  )
);

# Building command line
extra <- "vr = mutations";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep = ", "
  );
}

command <- base::paste0(
  "SomaticSignatures::motifMatrix(", extra, ")"
);
print(command)

# Building motif matrix
motifs <- base::eval(
  base::parse(
    text = command
  )
);

# Save results
if ("rds" %in% names(snakemake@output)) {
  base::saveRDS(
    object = motifs,
    file = base::as.character(x = snakemake@output[["rds"]])
  );
}

if ("tsv" %in% names(snakemake@output)) {
  motifs <- base::as.data.frame(motifs);
  utils::write.table(
    x = motifs,
    file = snakemake@output[["tsv"]],
    quote = FALSE,
    sep = "\t",
    row.names = TRUE
  );
}
