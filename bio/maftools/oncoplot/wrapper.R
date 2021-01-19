#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for maftools compute trinucleotideMatrix

# Many libraries are useless. This wrapper will be splitted in at least four
# of them.
base::library(package = "maftools", quietly = TRUE);

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

# Building trinucleotideMatrix computation command line
maf_obj <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
);

vc_nonSyn <- snakemake@params[["non_synonymous"]];
color_named_vector <- snakemake@params[["color_named_vector"]];
names(color_named_vector) <- vc_nonSyn;
oncoplot_extra <- "maf = maf_obj, color = color_named_vector";

if ("oncoplot_extra" %in% base::names(snakemake@params)) {
  oncoplot_extra <- base::paste(
    oncoplot_extra, snakemake@params[["oncoplot_extra"]], sep = ", "
  );
}
oncoplot_cmd_line <- base::paste0(
  "maftools::oncoplot(", oncoplot_extra, ")"
);

# Plotting signatures
if ("png" %in% base::names(snakemake@output)) {
  base::message(png_cmd_line);
  base::message(oncoplot_cmd_line);
  base::eval(base::parse(text = png_cmd_line));
  base::eval(base::parse(text = oncoplot_cmd_line));
  dev.off();
}
