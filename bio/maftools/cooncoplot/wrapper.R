#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for maftools co-oncoplot

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
if ("cooncoplot_extra" %in% base::names(snakemake@params)) {
  cooncoplot_extra <- base::paste(
    cooncoplot_extra,
    base::names(snakemake@params[["cooncoplot_extra"]]),
    sep=", "
  );
}
cooncoplot_cmd <- base::paste0(
  "maftools::mafCompare(", cooncoplot_extra, ")"
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
base::message(cooncoplot_cmd);

base::eval(base::parse(text=png_cmd_line));
base::eval(base::parse(text=cooncoplot_cmd));
dev.off();
