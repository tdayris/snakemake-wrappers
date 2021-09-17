#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for somatic interactions


base::library(package = "maftools", quietly = TRUE);
base::message("Libraries loaded");

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


# Building maftools command line
maf_obj <- base::readRDS(
  file = base::as.character(x = snakemake@input[["rds"]])
);
base::message("MAF loaded");

somatic_inter_extra <- "maf=maf_obj";
if ("somatic_inter_extra" %in% base::names(snakemake@params)) {
  somatic_inter_extra <- base::paste(
    somatic_inter_extra,
    base::as.character(x=snakemake@params[["extra"]]),
    sep=", "
  );
}

somtaic_interaction_cmd_line <- base::paste0(
  "maftools::somaticInteractions(", somatic_inter_extra, ")"
);
base::message(somtaic_interaction_cmd_line);

# Running command lines
base::eval(base::parse(text = png_cmd_line));
base::eval(base::parse(text = somtaic_interaction_cmd_line));
dev.off();
