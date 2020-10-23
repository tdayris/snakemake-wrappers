#!/usr/bin/R

# This script loads peaks into R with ChipSeeker

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"


# Load libraries
base::library(package = "ChIPseeker", quietly = TRUE);


# Load user parameters
txdb <- "hg38";
if ("txdb" %in% base::names(snakemake@params)) {
  txdb <- base::as.character(x = snakemake@params[["txdb"]]);
}

if (txdb == "hg19") {
  base::library(package = "TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE);
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene;
} else if (txdb == "hg38") {
  base::library(package = "TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE);
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
} else if (txdb == "mm10") {
  base::library(package = "TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE);
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene;
} else {
  base::library(package = "TxDb.Mmusculus.UCSC.mm9.knownGene", quietly = TRUE);
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene;
}

peak_files_path <- sapply(
  snakemake@input[["peaks"]],
  function(peak_path) base::as.character(x = peak_path)
);

rds_file_path <- base::as.character(
  x = snakemake@output[["rds"]]
);

# Load peak files
extra <- 'TxDb = txdb';
if ('get_promoters' %in% base::names(snakemake@params)) {
  cmd <- base::paste(
    extra,
    snakemake@params[["get_promoters"]],
    sep = ", "
  );
}

# Create object
promoters <- base::eval(
  base::parse(
    text = base::paste0(
      "ChIPseeker::getPromoters(", extra, ");"
    )
  )
);

tagMatrixList <- lapply(
  peak_files_path,
  ChIPseeker::getTagMatrix,
  windows=promoters
);


# Save RDS
base::saveRDS(
  object = tagMatrixList,
  file = rds_file_path
);
