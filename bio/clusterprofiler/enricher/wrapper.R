#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2022, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for ClusterProfiler enrichment analysis

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
        data <- clusterProfiler::read.gmt(gmtfile = data_path)
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

# Load gene information
genes <- read_input(data_path = snakemake@input[["gene"]])
weight <- stats::setNames(genes[, 2], genes[, 1])
genes <- base::unique(x = base::names(weight))

# Load set of enrichment terms
term2gene <- read_input(data_path = snakemake@input[["term2gene"]])
universe <- term2gene[, 2]

# Optional human readable terms
term2name <- NA
if ("term2name" %in% base::names(snakemake@input)) {
    term2name <- read_input(data_path = snakemake@input[["term2name"]])
}

# Build enricher function parameters
enrich_parameters <- "gene = genes, TERM2GENE = term2gene, universe = universe, TERM2NAME = term2name"
if ("enricher_extra" %in% base::names(snakemake@params)) {
    enrich_parameters <- base::paste(
        enrich_parameters,
        as.character(x = snakemake@params[["enricher_extra"]]),
        sep = ", "
    )
}

# Build and execute command line
enrich_command <- base::paste0(
    "clusterProfiler::enricher(",
    enrich_parameters,
    ")"
)
base::message(enrich_command)
enriched_terms <- base::eval(base::parse(text = enrich_command))

base::print(genes)
base::print(universe)
base::print(term2gene)
base::print(term2name)
base::print(enriched_terms)

# On user request, save enriched terms as RDS binary file.
if ("rds" %in% base::names(snakemake@output)) {
    base::saveRDS(
        object = enriched_terms, 
        file = base::as.character(x = snakemake@output[["rds"]])
    )
}


# On user request, save enriched terms as TSV formatted text file.
if ("tsv" %in% base::names(snakemake@output)) {
    utils::write.table(
        x = as.data.frame(x = enriched_terms),
        file = base::as.character(x = snakemake@output[["tsv"]]),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
    )
}


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()