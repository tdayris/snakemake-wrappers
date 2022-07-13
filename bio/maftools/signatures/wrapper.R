#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for maftools extractSignatures

base::library(package = "maftools", quietly = TRUE);
base::library(package = "NMF", quietly = TRUE);

# Building and running signature extration command line
maf.tnm <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
);
sig_extra <- "mat = maf.tnm";
if ("sig_extra" %in% base::names(snakemake@params)) {
  sig_extra <- base::paste(
    sig_extra, snakemake@params[["sig_extra"]], sep = ", "
  );
}
sig_cmd_line <- base::paste0("maftools::extractSignatures(", sig_extra, ")");
base::message(sig_cmd_line);
maf.sig <- base::eval(base::parse(text = sig_cmd_line));

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

# Building plotSignatures command line
plotsig_extra <- "nmfRes = maf.sig";
if ("plotsig_extra" %in% base::names(snakemake@params)) {
  plotsig_extra <- base::paste(
    plotsig_extra, snakemake@params[["plotsig_extra"]], sep = ", "
  );
}
plotsig_cmd_line <- base::paste0(
  "maftools::plotSignatures(", plotsig_extra, ")"
);

# Plotting signatures
if ("png" %in% base::names(snakemake@output)) {
  base::message(png_cmd_line);
  base::message(plotsig_cmd_line)
  base::eval(base::parse(text = png_cmd_line));
  base::eval(base::parse(text = plotsig_cmd_line));
  dev.off();
}

# Saving results
if ("rds" %in% base::names(snakemake@output)) {
  base::saveRDS(
    file = base::as.character(x = snakemake@output[["rds"]]),
    obj = maf.sig
  );
}
