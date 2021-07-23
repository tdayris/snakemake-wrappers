#!/usr/bin/env R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2021, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for ConsensusClusterPlus

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file <- base::file(snakemake@log[[1]], open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

mkdirs <- function(fp) {
    if(!base::file.exists(fp)) {
        mkdirs(base::dirname(fp));
        base::message("Creating ", fp);
        base::dir.create(fp);
    } else {
      base::message("Path already exists: ", fp);
    }
}

# Load input dataset
expr_mat <- utils::read.table(
  file = base::as.character(x=snakemake@input[["expr_mat"]]),
  header = TRUE,
  sep = "\t"
);
expr_mat <- base::as.matrix(expr_mat);
base::message(head(expr_mat));

# Build command line
extra_consensus <- "d=expr_mat, verbose=TRUE";
if ("csv" %in% base::names(snakemake@output)) {
  extra_consensus <- base::paste(
    extra_consensus,
    "writeTable=TRUE",
    sep=", "
  );
} else if ("res_dir" %in% base::names(snakemake@output)) {
  extra_consensus <- base::paste(
    extra_consensus,
    "writeTable=TRUE",
    sep=", "
  );

  mkdirs(
    fp=base::as.character(x=snakemake@output[["res_dir"]])
  );
  setwd(
    base::dirname(
      base::as.character(x=snakemake@output[["res_dir"]])
    )
  );
} else {
  extra_consensus <- base::paste(
    extra_consensus,
    "writeTable=FALSE",
    sep=", "
  );
}


ext = "png";
if ("plot" %in% base::names(snakemake@output)) {
  ext = base::unlist(base::strsplit(
    x=snakemake@output[["plot"]],
    split="\\."
  ))[2];
}
base::message("Saving images as: ", ext)
extra_consensus <- base::paste(
  extra_consensus,
  "plot = ext",
  sep=", "
);


if ("extra" %in% base::names(snakemake@params)) {
  extra_consensus <- base::paste(
    extra_consensus,
    base::as.character(x=snakemake@params[["extra"]]),
    sep=", "
  );
}

ConsensusClusterPlus_cmd <- base::paste0(
  "ConsensusClusterPlus::ConsensusClusterPlus(",
  extra_consensus,
  ")"
);
base::message(ConsensusClusterPlus_cmd);

base::eval(
  base::parse(
    text = ConsensusClusterPlus_cmd
  )
);


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
