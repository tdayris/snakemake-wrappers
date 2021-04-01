#!/bin/R

# Loading library
base::library("tximport");   # Perform actual count importation in R
base::library("readr");      # Read faster!
base::library("jsonlite");   # Importing inferential replicates

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

# Cast input paths as character to avoid errors
samples_paths <- sapply(               # Sequentially apply
  snakemake@input[["quant"]],          # ... to all quantification paths
  function(quant) as.character(quant)  # ... a cast as character
);

# Collapse path into a character vector
samples_paths <- base::paste0(samples_paths, collapse = '", "');
if ("sample_names" %in% base::names(snakemake@params)) {
  names(samples_paths) <- snakemake@params[["names"]];
}

# Building function arguments
extra <- base::paste0('files = c("', samples_paths, '")');

# Check if user provided optional transcript to gene table
if ("tx_to_gene" %in% names(snakemake@input)) {
  tx2gene <- readr::read_tsv(snakemake@input[["tx_to_gene"]]);
  extra <- base::paste(
    extra,                 # Foreward existing arguments
    ", tx2gene = ",        # Argument name
    "tx2gene"              # Add tx2gene to parameters
  );
}

# Add user defined arguments
if ("extra" %in% names(snakemake@params)) {
  if (snakemake@params[["extra"]] != "") {
    extra <- base::paste(
      extra,                       # Foreward existing parameters
      snakemake@params[["extra"]], # Add user parameters
      sep = ", "                   # Field separator
    );
  }
}


print(extra);
# Perform tximport work
txi <- base::eval(                        # Evaluate the following
  base::parse(                            # ... parsed expression
    text = base::paste0(
      "tximport::tximport(", extra, ");"  # ... of tximport and its arguments
    )
  )
);

# Save results
base::saveRDS(                       # Save R object
  object = txi,                      # The txi object
  file = snakemake@output[["txi"]]   # Output path is provided by Snakemake
);

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
