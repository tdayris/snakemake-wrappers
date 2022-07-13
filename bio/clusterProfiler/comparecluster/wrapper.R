#!/usr/bin/env R

# This script takes a list of geneList objects and performs
# an analysis on each one based on statistical formula

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2022, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");


readGeneList <- function(path) {
    genelist <- utils::read.table(
        file=base::as.character(x=path), 
        header=TRUE, 
        sep="\t"
    );
    return(genelist)
}

readGeneLists <- function(paths, conditions) {
    # Load all TSV files into one dataframe
    genelists <- base::lapply(
        paths, function(p) readGeneList(p)
    )
    # Add condition and return the dataframe
    genelists$Condition <- condition
    return(genelists)
}

conditions <- sapply(
    snakemake@params[["condition"]]),
    function(c) base::as.character(x=c)
);

geneLists <- readGeneLists(
    paths=snakemake@input[["gene_list"]],
    conditions=conditions
);

formula <- stats::as.formula(x="~Condition");
paramfunction <- base::as.character(
    x=snakemake@params[["fun"]]
);

extra <- "formula, data=geneLists, fun=paramfunction";
if ("extra" %in% base::names(snakemake@params)) {
    extra <- base::paste(
        extra,
        base::as.character(x=snakemake@params[["extra"]]),
        sep=", "
    );
}

command <- base::paste0(
    "clusterProfiler::compareCluster(",
    extra,
    ")"
);