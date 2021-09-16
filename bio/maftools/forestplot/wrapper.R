#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper that compares two maf cohorts

base::library(package="maftools", quietly=TRUE);


maf_paths <- sapply(
  mafrds, function(p) base::as.character(x=p)
);

maf_names <- sapply(
  mafnames, function(n) base::as.character(x=n)
);

maf_objs <- sapply(
  maf, function(m) base::readRDS(m)
);
base::message("MAF loaded");

compare_extra <- "m1=maf_obj[1], m2=maf_obj[2], m1Name=maf_names[1], m2Name=maf_names[2]";
if ("compare_extra" %in% base::names(snakemake@params)) {
  compare_extra <- base::paste(
    compare_extra,
    base::names(snakemake@params[["compare_extra"]]),
    sep=", "
  );
}
compare_cmd <- base::paste0(
  "maftools::mafCompare(", compare_extra, ")"
);
base::message(compare_cmd);

compared <- base::eval(base::parse(text=compare_cmd));

if ("tsv" %in% base::names(snakemake@output)) {
  out_tsv <- base::as.character(x=snakemake@output[["tsv"]]);

  utils::write.table(
    base::as.data.frame(compared),
    sep="\t",
    row.names=FALSE,
    col.names=TRUE
  );
}

if ("png" %in% base::names(snakemake@output)) {
  out_png <- base::as.character(x=snakemake@output[["png"]]);

  forest_extra <- "mafCompareRes=compared";
  if ("forest_extra" %in% base::names(snakemake@params)) {
    forest_extra <- base::paste(
      forest_extra,
      base::as.character(x=snakemake@params[["forest_extra"]]),
      sep=", "
    );
  }

  forest_cmd <- base::paste(
    "maftools::forestPlot(", forest_extra, ")"
  );

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
  base::message(forest_cmd);

  base::eval(base::parse(text=png_cmd_line));
  base::eval(base::parse(text=forest_cmd));
  dev.off();
}
