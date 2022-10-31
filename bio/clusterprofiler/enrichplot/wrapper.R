#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2022, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Snakemake wrapper for ClusterProfiler enrichment analysis

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file <- file(snakemake@log[[1]], open = "wt")
base::sink(log.file)
base::sink(log.file, type = "message")

# Load packages
base::library(package = "ggplot2", character.only = TRUE)
base::library(package = "clusterProfiler", character.only = TRUE)
base::library(package = "Cairo", character.only = TRUE)
base::library(package = "UpSetR", character.only = TRUE)
base::library(package = "enrichplot", character.only = TRUE)



# Read snakemake input, and guess wether it's a RDS, a TSV, or a GMT file
read_input <- function(data_path) {
    base::message(base::paste("Loading", data_path))

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

# This function performs graphs.
plot_enrichment <- function(plot_function, plot_base_extra, output_plot_key, param_plot_key) {
    # Acquire grDevices::png parameters
    out_png <- base::as.character(x = snakemake@output[[output_plot_key]])
    base::message("Saving plot to ", out_png)
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
    base::print(base::eval(base::parse(text = plot_command)))
    grDevices::dev.off()
}

# Load enrichment/gsea input data
enrichment <- read_input(
    data_path = base::as.character(x = snakemake@input[["rds"]])
)
base::print(head(enrichment))

# In order to color heatplot, provide gene weights
genes <- NULL
weights <- NULL
if ("genes" %in% base::names(snakemake@input)) {
    genes <- read_input(
        data_path = base::as.character(x = snakemake@input[["genes"]])
    )
    weights <- stats::setNames(genes[, 2], genes[, 1])
}

# In order to plot emapplot, one need to simplify terms with
# a GO semantic similarity object
go_sem_sim_object <- NULL
if ("go_sem_sim" %in% base::names(snakemake@input)) {
    go_sem_sim_object <- read_input(
        data_path = base::as.character(x = snakemake@input[["go_sem_sim"]])
    )
}

# On user request, save barplot
if ("barplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "graphics::barplot",
        plot_base_extra = "enrichment",
        output_plot_key = "barplot",
        param_plot_key = "barplot_extra"
    )
}

# On user request, save dotplot
if ("dotplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "enrichplot::dotplot",
        plot_base_extra = "object = enrichment",
        output_plot_key = "dotplot",
        param_plot_key = "dotplot_extra"
    )
}

# On user request, save cnetplot
if ("cnetplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "enrichplot::cnetplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "cnetplot",
        param_plot_key = "cnetplot_extra"
    )
}

# On user request, save heatmap
if ("heatplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "enrichplot::heatplot",
        plot_base_extra = "x = enrichment, foldChange = weights",
        output_plot_key = "heatplot",
        param_plot_key = "heatplot_extra"
    )
}

# On user request, save emapplot
if ("emapplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "enrichplot::emapplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "emapplot",
        param_plot_key = "emapplot_extra"
    )
}

# On user request, save upsetplot
if ("upsetplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "enrichplot::upsetplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "upsetplot",
        param_plot_key = "upsetplot_extra"
    )
}


# On user request, save pmcplot
if ("pmcplot" %in% base::names(snakemake@output)) {
    terms <- enrichment$Description
    plot_enrichment(
        plot_function = "enrichplot::pmcplot",
        plot_base_extra = "query = terms",
        output_plot_key = "pmcplot",
        param_plot_key = "pmcplot_extra"
    )
}

# On user request, save a ridgeplot
if ("ridgeplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "enrichplot::ridgeplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "ridgeplot",
        param_plot_key = "ridgeplot_extra"
    )
}

# On user request, save a gsea plot
if ("gseaplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "enrichplot::gseaplot2",
        plot_base_extra = "x = enrichment",
        output_plot_key = "gseaplot",
        param_plot_key = "gseaplot_extra"
    )
}

if ("gsearank" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "enrichplot::gsearank",
        plot_base_extra = "x = enrichment",
        output_plot_key = "gsearank",
        param_plot_key = "gsearank_extra"
    )
}

if ("emapplot" %in% base::names(snakemake@output)) {
    pairwise_termism_params <- extra_parameters(
        parameters = "x = enrichment, semData = go_sem_sim_object",
        param_key = "pairwise_termism"
    )
    pairwise_termism_command <- paste0(
        "enrichplot::pairwise_termism(",
        pairwise_termism_params,
        ")"
    )
    base::message(pairwise_termism_command)
    enrichment <- base::eval(base::parse(text = pairwise_termism_command))

    plot_enrichment(
        plot_function = "enrichplot::emapplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "emapplot",
        param_plot_key = "emapplot_extra"
    )
}


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()