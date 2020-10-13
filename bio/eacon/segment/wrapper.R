#!/bin/R

base::library(package = "EaCoN", quietly = TRUE);

input_rds <- base::as.character(
  x = snakemake@input[["rds"]]
);

extra <- "force=TRUE";
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::as.character(x = snakemake@params[["extra"]]);
}

command <- base::paste0(
  "EaCoN::Segment.ff(RDS.file = input_rds, ",
  extra,
  ")"
);
base::message(command);

base::eval(
  base::parse(
    text=command
  )
);
