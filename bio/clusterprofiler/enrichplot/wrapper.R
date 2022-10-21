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
base::library(package = "clusterProfiler", character.only = TRUE)
base::library(package = "Cairo", character.only = TRUE)
base::library(package = "UpSetR", character.only = TRUE)
base::library(package = "enrichplot", character.only = TRUE)


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
    base::eval(base::parse(text = plot_command))
    grDevices::dev.off()
}


# Load enrichment/gsea input data
enrichment <- base::readRDS(
    file = base::as.character(x = snakemake@input[["rds"]])
)
base::print(head(enrichment))



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
        plot_function = "dotplot",
        plot_base_extra = "object = enrichment",
        output_plot_key = "dotplot",
        param_plot_key = "dotplot_extra"
    )
}

# On user request, save cnetplot
if ("cnetplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "cnetplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "cnetplot",
        param_plot_key = "cnetplot_extra"
    )
}

# On user request, save heatmap
if ("heatplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "heatplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "heatplot",
        param_plot_key = "heatplot_extra"
    )
}

# On user request, save emapplot
if ("emapplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "emapplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "emapplot",
        param_plot_key = "emapplot_extra"
    )
}

# On user request, save upsetplot
if ("upsetplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "upsetplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "upsetplot",
        param_plot_key = "upsetplot_extra"
    )
}


# On user request, save pmcplot
if ("pmcplot" %in% base::names(snakemake@output)) {
    terms <- enrichment$Description
    plot_enrichment(
        plot_function = "pmcplot",
        plot_base_extra = "query = terms",
        output_plot_key = "pmcplot",
        param_plot_key = "pmcplot_extra"
    )
}

# On user request, save a ridgeplot
if ("ridgeplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "ridgeplot",
        plot_base_extra = "x = enrichment",
        output_plot_key = "ridgeplot",
        param_plot_key = "ridgeplot_extra"
    )
}

# On user request, save a gsea plot
if ("gseaplot" %in% base::names(snakemake@output)) {
    plot_enrichment(
        plot_function = "gseaplot2",
        plot_base_extra = "x = enrichment",
        output_plot_key = "gseaplot"
    )
}


# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()