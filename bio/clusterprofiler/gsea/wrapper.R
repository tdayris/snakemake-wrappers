#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2022, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for ClusterProfiler gsea analysis

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]], open = "wt")
base::sink(log.file)
base::sink(log.file, type = "message")

# Load packages
base::library(package = "clusterProfiler", character.only = TRUE)


# Read snakemake input, and guess wether it's a RDS, a TSV, or a GMT file
read_input <- function(data_path) {
    data_path <- base::as.character(x = data_path)

    data <- NULL
    if (base::endsWith(x = data_path, suffix = ".RDS")) {
        # Then the binary RDS file is loaded
        data <- readRDS(file = data_path)
    } else if (base::endsWith(x = data_path, suffix = ".gmt")) {
        # Then the file is a GMT file (e.g. from MSigDB)
        data <- clusterProfiler::read.gmt(gmtfile = gmt_path)
    } else if (base::endsWith(x = data_path, suffix = ".csv")) {
        # Then the file is a CSV
        data <- utils::read.table(
            file = data_path,
            header = FALSE,
            sep = ",",
            stringsAsFactors = FALSE
        )
    } else {
        # Then it is a TSV text file, since it's the only
        # other expected format.
        data <- utils::read.table(
            file = data_path,
            header = FALSE,
            sep = "\t",
            stringsAsFactors = FALSE
        )
    }

    base::return(data)
}

# This function creates ClusterProfiler gene_list object
build_gene_list <- function(gene_data_frame) {
    # ClusterProfiler requires its input dataset (aka
    # `gene_list`) to be named vectors. This function
    # turn dataframes into named vectors.
    genes_names_vector <- stats::setNames(gene_data_frame[, 2], gene_data_frame[, 1])
    base::print(genes_names_vector)
    base::return(genes_names_vector)

}

# This function gathers optional parameters
extra_parameters <- function(parameters, param_key) {
    # if optional parameters are provided by use
    # in snakemake@params, then add them to the
    # already existing parameters
    if (param_key %in% base::names(snakemake@params)) {
        parameters <- base::paste(
            parameters,
            as.character(x = snakemake@params[[param_key]]),
            sep = ", "
        )
    }

    base::return(parameters)
}

# Load gene information
genes <- read_input(data_path = snakemake@input[["gene"]])
weight <- stats::setNames(genes[, 2], genes[, 1])
weight <- base::sort(x = weight, decreasing = TRUE)

# Load set of gseament terms
term2gene <- read_input(data_path = snakemake@input[["term2gene"]])

# Build gsea function parameters
gsea_parameters <- "geneList = weight, TERM2GENE = term2gene"
if ("gsea_extra" %in% base::names(snakemake@params)) {
    gsea_parameters <- base::paste(
        gsea_parameters,
        as.character(x = snakemake@params[["gsea_extra"]]),
        sep = ", "
    )
}

# Add optional human-readable term names
term2name <- NA
if ("term2name" %in% base::names(snakemake@input)) {
    term2name <- read_input(data_path = snakemake@input[["term2name"]])
    gsea_parameters <- base::paste(
        gsea_parameters,
        "TERM2NAME = term2name",
        sep = ", "
    )
}

# Build and execute command line
gsea_command <- base::paste0(
    "clusterProfiler::GSEA(",
    gsea_parameters,
    ")"
)
base::message(gsea_command)
gsea_terms <- base::eval(base::parse(text = gsea_command))


# On user request, save gsea terms as RDS binary file.
if ("rds" %in% base::names(snakemake@output)) {
    base::saveRDS(
        object = gsea_terms,
        file = base::as.character(x = snakemake@output[["rds"]])
    )
}


# On user request, save gsea terms as TSV formatted text file.
if ("tsv" %in% base::names(snakemake@output)) {
    utils::write.table(
        x = as.data.frame(x = gsea_terms),
        file = base::as.character(x = snakemake@output[["tsv"]]),
        sep = "\t",
    )
}


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()