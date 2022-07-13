#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for maftools compute trinucleotideMatrix

# Many libraries are useless. This wrapper will be splitted in at least four
# of them.
base::library(package = "maftools", quietly = TRUE);
base::library(package = "BSgenome.Mmusculus.UCSC.mm9", quietly = TRUE);
base::library(package = "NMF", quietly = TRUE);

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
trinuc_mat_extra <- "maf = maf_obj";
if ("trinuc_mat_extra" %in% base::names(snakemake@params)) {
  trinuc_mat_extra <- base::paste(
    trinuc_mat_extra,
    "ref_genome = 'BSgenome.Mmusculus.UCSC.mm9'",
    snakemake@params[["trinuc_mat_extra"]],
    sep = ", "
  );
}

trinucleotide_cmd_line <- base::paste0(
  "maftools::trinucleotideMatrix(", trinuc_mat_extra, ")"
);

# Building cophrenic correlation coefficient plot command line
# The results of trinucleotideMatrix will be stored in a
# variable called: maf.tnm ; don't be surprised on this poping variable name.
estimate_extra <- "mat = maf.tnm";
if ("estimate_extra" %in% base::names(snakemake@params)) {
  estimate_extra <- base::paste(
    estimate_extra, snakemake@params[["estimate_extra"]], sep = ", "
  );
}
estimate_cmd_line <- base::paste0(
  "maftools::estimateSignatures(",
  estimate_extra,
  ")"
);

# Running command lines on user's request
base::message(trinucleotide_cmd_line);
maf.tnm <- base::eval(base::parse(text = trinucleotide_cmd_line));

if ("tsv" %in% base::names(snakemake@output)) {
  utils::write.table(
    x = maf.tnm,
    file = base::as.character(x = snakemake@output[["tsv"]]),
    sep = "\t",
    quote = FALSE
  );
}

if ("png" %in% base::names(snakemake@output)) {
  base::message(png_cmd_line);
  base::eval(base::parse(text = png_cmd_line));

  base::message(estimate_cmd_line);
  base::eval(base::parse(text = estimate_cmd_line));
  dev.off();
}

if ("rds" %in% base::names(snakemake@output)) {
  base::saveRDS(
    file = base::as.character(x = snakemake@output[["rds"]]),
    obj = maf.tnm
  );
}
