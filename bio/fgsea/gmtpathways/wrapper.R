#!/usr/bin/R

# This script takes a GMT file and loads it into R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2021, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

base::library(package = "fgsea");

gmt_file <- as.character(x = snakemake@input[["gmt"]]);

gmt_data <- fgsea::gmtPathways(gmt.file = gmt_file);

base::saveRDS(
  object = gmt_data,
  file = base::as.character(snakemake@output[["rds"]])
);

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
