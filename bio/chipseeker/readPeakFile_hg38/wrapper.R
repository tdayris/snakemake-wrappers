#!/usr/bin/R

# This script loads peaks into R with ChipSeeker

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2021, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Load libraries
base::library(package = "ChIPseeker", quietly = TRUE);
base::library(package = "TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE);


# Cast as character since '/' causes issues in R
bed_file <- base::as.character(x=snakemake@input[["bed"]]);
extra_readpeaks <- "bed_file, as = 'GRanges'";
if ("extra_readpeaks" %in% base::names(snakemake@params)) {
  extra_readpeaks <- base::paste(
    extra_readpeaks,
    base::as.character(x=snakemake@params[["extra_readpeaks"]]),
    sep=","
  );
}

extra_promoters <- "TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene";
if ("extra_promoters" %in% base::names(snakemake@params)) {
  extra_promoters <- base::paste(
    extra_promoters,
    base::as.character(x=snakemake@params[["extra_promoters"]]),
    sep=","
  );
}

extra_tagmatrix <- "peak=granges, windows=promoters";
if ("extra_tagmatrix" %in% base::names(snakemake@params)) {
  extra_tagmatrix <- base::paste(
    extra_tagmatrix,
    base::as.character(x=snakemake@params[["extra_tagmatrix"]]),
    sep=","
  );
}


# Build GRanges, getPromoters and TagMatrix command lines
granges_command <- base::paste0("ChIPseeker::ReadPeakFile(", extra, ");");
base::messages(granges_command);
get_promoters_command <- base::paste0("ChIPseeker::getPromoters(", extra, ")");
base::messages(get_promoters_command);
tagmatrix_command <- base::paste0("ChIPseeker::getTagMatrix(", extra, ");");
base::messages(tagmatrix_command);

# Build GRanges and TagMatrix objects
granges <- base::eval(base::parse(text = granges_command));
promoters <- base::eval(base::parse(text = get_promoters_command));
tagmatrix <- base::eval(base::parse(text = tagmatrix_command));


# Keeping results
if ("grange" %in% base::names(snakemake@output)) {
  base::saveRDS(
    object=granges,
    file=base::as.character(x=snakemake@output[["grange"]])
  );
}

if ("tagmatrix" %in% base::names(snakemake@output)) {
  base::saveRDS(
    object=tagmatrix,
    file=base::as.character(x=snakemake@output[["tagmatrix"]])
  );
}
