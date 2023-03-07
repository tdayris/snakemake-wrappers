#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log.file)
base::sink(log.file, type = "message")

# Load Libraries
base::library(package = "dplyr", character.only = TRUE)
base::library(package = "tidyr", character.only = TRUE)
base::library(package = "tibble", character.only = TRUE)
base::library(package = "readr", character.only = TRUE)
base::library(package = "ggplot2", character.only = TRUE)
base::library(package = "immunedeconv", character.only = TRUE)

base::message("Libraries loaded")

# Load dataset
gene_col <- "GENE";
if ("gene_col" %in% base::names(snakemake@params)) {
  gene_col <- base::as.character(x = snakemake@params[["gene_col"]])
}

tpm <- utils::read.table(
  file = base::as.character(x = snakemake@input[["expr_mat"]]),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
rownames(tpm) <- tpm[, gene_col]
tpm[, gene_col] <- NULL
print(tpm %>% head)

method <- "quantiseq"
if ("method" %in% base::names(snakemake@params)) {
  method <- base::as.character(x = snakemake@params[["method"]])
}

extra <- "method = method"
if ("extra" %in%base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    base::as.character(x = snakemake@params[["extra"]]),
    sep = ", "
  )
}

if (method == "cibersort" || method == "cibersort_abs") {
  if ("cibersort_binary" %in% base::names(snakemake@input)) {
    immunedeconv::set_cibersort_binary(
      base::as.character(x = snakemake@input[["cibersort_binary"]])
    )
  } else {
    base::message("## WARNING: No cibersort binary provided ##")
  }

  if ("cibersort_mat" %in% base::names(snakemake@input)) {
    immunedeconv::set_cibersort_mat(
      base::as.character(x = snakemake@input[["cibersort_mat"]])
    )
  } else {
    base::message("## WARNING: No cibersort matrix provided ##")
  }
}

cmd <- base::paste0(
  "immunedeconv::deconvolute(",
  "gene_expression = tpm, ",
  extra,
  ")"
)
base::message("Datasets and configuration loaded")
print(cmd)

# Deconvolution
res_deconv <- base::eval(base::parse(text = cmd))
colors <- rainbow(
  n = length(unique(unlist(res_deconv["cell_type"]))),
  s = 0.75,
  v = 1,
)
print(res_deconv %>% head)
base::message("Deconvolution performed")

# Save results
if ("rds" %in% base::names(snakemake@output)) {
  base::saveRDS(
    obj = res_deconv,
    file = snakemake@output[["rds"]]
  )
  base::message("RDS object saved as ", snakemake@output[["rds"]])
}

if ("tsv" %in% base::names(snakemake@output)) {
  utils::write.table(
    x = res_deconv,
    file = snakemake@output[["tsv"]],
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  base::message("TSV table saved as ", snakemake@output[["tsv"]])
}

# Plot graphs
if ("histogram" %in% base::names(snakemake@output)) {
  png(
    filename = snakemake@output[["histogram"]],
    width = dotx,
    height = doty,
    units = "px",
    type = "cairo"
  )

  print(res_deconv %>%
    gather(sample, fraction, -cell_type) %>%
    ggplot(aes(x = sample, y = fraction, fill = cell_type)) +
      geom_bar(stat='identity') +
      coord_flip() +
      scale_fill_manual(values = colors) +
      scale_x_discrete(
        limits = rev(levels(res_deconv))
      )
    )

  dev.off()
  base::message("Histogram saved as ", snakemake@output[["histogram"]])
}


if ("dotplot" %in% base::names(snakemake@output)) {
  png(
    filename = snakemake@output[["dotplot"]],
    width = dotx,
    height = doty * 4,
    units = "px",
    type = "cairo"
  )

  print(res_deconv %>%
    gather(sample, score, -cell_type) %>%
    ggplot(aes(x = sample, y = score)) +
      geom_point(size = 4) +
      facet_wrap(~cell_type, scales = "free_x", ncol = 3) +
      coord_flip() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
    )

  dev.off()
  base::message("Dotplots saved as ", snakemake@output[["dotplot"]])
}

if ("plotdir" %in% base::names(snakemake@output)) {
  plotdir = base::as.character(x = snakemake@output[["plotdir"]])
  if (! dir.exists(snakemake@output[["plotdir"]])) {
    dir.create(plotdir)
  }
  for (celltype in base::unlist(res_deconv["cell_type"])) {
    png_name <- base::paste(
      base::gsub(" ", "_", celltype), "dotplot", "png", sep = "."
    )
    png_path <- base::file.path(plotdir, png_name);
    base::message("Saving ", celltype, " dotplot, as: ", png_path)

    png(
      filename = png_path,
      width = dotx,
      height = doty,
      units = "px",
      type = "cairo"
    )

    print(res_deconv %>%
      filter(cell_type == celltype) %>%
      gather(sample, score, -cell_type) %>%
      ggplot(aes(x = sample, y = score)) +
        geom_point(size = 4) +
        coord_flip() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    )

    dev.off()
  }
}

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()