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
base::library(package = "clusterProfiler", character.only = TRUE);
base::library(package = "Cairo", character.only = TRUE);
base::library(package = "enrichplot", character.only = TRUE);


# Read snakemake input, and guess wether it's a RDS, a TSV, or a GMT file
read_input <- function(data_path) {
    data_path <- base::as.character(x = data_path)

    data <- NULL
    if (base::endsWith(x = data_path, suffix = ".RDS") {
        # Then the binary RDS file is loaded
        data <- readRDS(file = data_path)
    } else if (base::endsWith(x = data_path, suffix = ".gmt")) {
        # The the file is a GMT file (e.g. from MSigDB)
        data <- clusterProfiler::read.gmt(gmtfile = gmt_path)
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
    gene_list <- gene_data_frame[, 1]
    gene_list <- base::as.numeric(gene_list)

    # Note ENTREZ-genes id should be all numeric in case of
    # intersection with databases like org.XX.eg.db or
    # DOSE, KEGG, NCG, ...
    # If any of the identifier is not full numeric,
    # then genes identifiers are not ENTREZ-gene id,
    # and should be left as they are.
    genes_id <- gene_data_frame[, 0]
    genes_id <- base::tryCatch(
        expr = base::as.numeric(x = gene_data_frame[, 0]),
        error = function(e) { base::return(gene_data_frame[, 0]) },
        warning = function(e) { base::return(gene_data_frame[, 0]) },
    )

    base::return(gene_list)
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

# This function performs optional graphs on user demand.
plot_enrichment <- function(plot_function, plot_base_extra, output_plot_key, param_plot_key) {
    # Acquire grDevices::png parameters
    out_png <- base::as.character(x = snakemake@output[[output_plot_key]])
    png_params <- extra_parameters(
        parameters = "filename = out_png",
        param_key = "png_extra"
    )

    # Acquire plot parameters
    plot_params <- extra_parameters(
        parameters = plot_base_extra,
        param_key = param_plot_key
    )

    # Build command lines
    png_command <- base::paste0("grDevices::png(", png_params, ")")
    base::message(png_command)

    plot_command <- base::paste0(plot_function, "(", plot_params, ")")
    base::message(plot_command)

    base::eval(base::parse(text = png_command))
    base::eval(base::parse(text = plot_command))
    grDevices::dev.off()
}

# Load gene information
weight <- build_gene_list(
    gene_data_frame = read_input(data_path = snakemake@input[["gene"]])
)
genes <- base::names(weight)

# Load set of enrichment terms
term2gene <- read_input(data_path = snakemake@input[["term2gene"]])
universe <- term2gene[, 1]

# Build enricher function parameters
enrich_parameters <- extra_parameters(
    paramters = "gene = genes, TERM2GENE = term2gene, universe = universe", 
    param_key = "enrich_extra"
)

# Add optional human-readable term names
term2name <- NULL
if ("term2name" %in% base::names(snakemake@input)) {
    term2name <- read_input(data_path = snakemake@input[["term2name"]])
    enrich_parameters <- base::paste(
        enrich_parameters,
        "TERM2NAME = term2name",
        sep = ", "
    )
}

# Build and execute command line
enrich_command <- base::paste0("clusterProfiler::enricher(", extra, ")")
base::message(enrich_command)
enriched_terms <- base::eval(base::parse(text=enrich_command))


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
    )
}



# On user request, save barplot
if ("barplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "barplot", 
        plot_base_extra = "enriched_terms", 
        output_plot_key = "barplot", 
        param_plot_key = "barplot_extra"
    )
}

# On user request, save dotplot
if ("dotplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "dotplot", 
        plot_base_extra = "object = enriched_terms", 
        output_plot_key = "dotplot", 
        param_plot_key = "dotplot_extra"
    )
}

# On user request, save cnetplot
if ("cnetplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "cnetplot", 
        plot_base_extra = "x = enriched_terms", 
        output_plot_key = "cnetplot", 
        param_plot_key = "cnetplot_extra"
    )
}

# On user request, save heatmap
if ("heatplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "heatplot",
        plot_base_extra = "x = enriched_terms, foldChange = weights",
        output_plot_key = "heatplot", 
        param_plot_key = "heatplot_extra"
    )
}

# On user request, save emapplot
if ("emapplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "emapplot",
        plot_base_extra = "x = enriched_terms",
        output_plot_key = "emapplot", 
        param_plot_key = "emapplot_extra"
    )
}

# On user request, save upsetplot
if ("upsetplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "upsetplot",
        plot_base_extra = "x = enriched_terms",
        output_plot_key = "upsetplot", 
        param_plot_key = "upsetplot_extra"
    )
}


# On user request, save pmcplot
if ("pmcplot" %in% base::names(snakemake@output)) {
    terms <- enriched_terms$Description
    plot_enrichment(
        plot_function = "pmcplot",
        plot_base_extra = "query = terms",
        output_plot_key = "pmcplot", 
        param_plot_key = "pmcplot_extra"
    )
}


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()