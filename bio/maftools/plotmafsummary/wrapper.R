#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for maftools plotmafSummary

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
base::message(png_cmd_line);

# Building maftools::plotmafSummary command
maf_obj <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
);

maftools_extra <- "maf = maf_obj";
if ("maftools_extra" %in% base::names(snakemake@params)) {
  maftools_extra <- paste(
    maftools_extra, snakemake@params[["maftools_extra"]], sep = ", "
  );
}
maf_cmd_line <- base::paste0("maftools::plotmafSummary(", maftools_extra, ")");
base::message(maf_cmd_line);

# Running both command lines
base::eval(base::parse(text = png_cmd_line));
base::eval(base::parse(text = maf_cmd_line));
dev.off();
