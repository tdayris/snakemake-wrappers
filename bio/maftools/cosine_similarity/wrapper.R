#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for both maftools compareSignatures,
# and pheatmap cosine similarity

base::library(package = "maftools", quietly = TRUE);
base::library(package = 'pheatmap', quietly = TRUE);
base::library(package = "NMF", quietly = TRUE);

# Building and running signature comparison
maf.sig <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
)
comp_sig_extra <- "nmfRes = maf.sig";
if ("comp_sig_extra" %in% base::names(snakemake@params)) {
  comp_sig_extra <- base::paste(
    comp_sig_extra, snakemake@params[["comp_sig_extra"]], sep = ", "
  );
}
sig_cmd_line <- base::paste0(
  "maftools::compareSignatures(",
  comp_sig_extra,
  ")"
);
maf.comp <- base::eval(base::parse(text = sig_cmd_line));
head(maf.comp)

# Building graphics environment command
png_extra <- base::paste0(
  "filename", "='", base::as.character(x = snakemake@output[["png"]]), "'"
);
if ("png_extra" %in% base::names(snakemake@params)) {
  png_extra <- base::paste(
    png_extra, snakemake@params[["png_extra"]], sep = ", "
  );
}
png_cmd_line <- base::paste0("grDevices::png(", png_extra, ")");
base::message(png_cmd_line);
base::eval(base::parse(text = png_cmd_line));

# Plotting cosine similarities
pheatmap::pheatmap(
  mat = maf.comp$cosine_similarities,
  cluster_rows = FALSE,
  main = "Cosine similarity against validated signatures"
);
dev.off();

# Saving results
if ("tsv" %in% base::names(snakemake@output)) {
  utils::write.table(
    x = maf.comp,
    file = base::as.character(x = snakemake@output[["tsv"]]),
    sep = "\t",
    quote = FALSE
  );
}

if ("rds" %in% base::names(snakemake@output)) {
  base::saveRDS(
    file = base::as.character(x = snakemake@output[["rds"]]),
    obj = maf.comp
  );
}
