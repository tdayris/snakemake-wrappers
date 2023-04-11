# This script takes a deseq2 dataset object and performs
# a default DESeq2 analysis

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file <- base::file(description = snakemake@log[[1]], open = "wt")
base::sink(file = log.file)
base::sink(file = log.file, type = "message")

# Differential Gene expression
base::library(package = "BiocParallel", character.only = TRUE)
base::library(package = "SummarizedExperiment", character.only = TRUE)
base::library(package = "DESeq2", character.only = TRUE)

# Setting up multithreading if required
parallel <- FALSE
if (snakemake@threads > 1) {
    BiocParallel::register(
      BPPARAM=BiocParallel::MulticoreParam(snakemake@threads)
    )
    parallel <- TRUE
}

# Load DESeq2 dataset
dds_path <- base::as.character(x = snakemake@input[["dds"]])
dds <- base::readRDS(file = dds_path)

# Ensuring stats model is the one expected by user.
if ("factor" %in% base::names(x = snakemake@params)) {
  # Then user provided a factor and relevels are possible.
  if ("reference_level" %in% names(x = snakemake@params)) {
    reference <- base::as.character(
      x = snakemake@params[["reference_level"]]
    )
    dds[[factor]] <- stats::relevel(dds[[factor]], ref = reference)
    base::message("Reference level is ", reference)
  }
}

# Build extra parameters for DESeq2
extra_deseq2 <- "object = dds"
if ("extra" %in% snakemake@params) {
  extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["extra"]])
  )
}
base::message("Libraries and dataset loaded")

deseq2_cmd <- base::paste0(
  "DESeq2::DESeq(", extra_deseq2, ")"
)
base::message("DESeq2 command line:")
base::message(deseq2_cmd)

# Running DESeq2::DESeq for wald test result
wald <- base::eval(base::parse(text = deseq2_cmd))


# Save results RDS on user request
if ("wald_rds" %in% base::names(x = snakemake@output)) {
  output_rds <- base::as.character(x = snakemake@output[["wald_rds"]])
  base::saveRDS(obj = wald, file = output_rds)
  base::message("Wald test over, RDS saved")
)

# Saving normalized counts on demand
table <- counts(wald)
if ("normalized_counts_table" %in% base::names(snakemake@output)) {
  output_table <- base::as.character(x=snakemake@output[["normalized_counts_table"]])
  utils::write.table(x = table, file = output_table, sep = "\t", quote = FALSE)
  base::message("Normalized counts saved as TSV")
}

if ("normalized_counts_rds" %in% base::names(snakemake@output)) {
  output_rds <- base::as.character(
    x = snakemake@output[["normalized_counts_rds"]]
  )
  base::saveRDS(
    obj = table,
    file = output_rds
  )
}

# On user request: saving all results as TSV in a directory.
# User can later access the directory content with
# a checkpoint rule.
if ("deseq2_result_dir" %in% base::names(snakemake@output)) {
  wald_results_names <- DESeq2::resultsNames(object = wald)

  # Recovreing extra parameters for TSV tables
  # The variable `result_name` is built below in `for` loop.
  extra_results <- "object = wald, name = result_name"
  if ("extra_results" %in% base::names(snakemake@params)) {
    extra_results <- base::paste(
      ", ",
      base::as.character(x = snakemake@params[["extra_results"]])
    )
  }

  # DESeq2 result dir will contain all results available in the Wald object
  output_prefix <- snakemake@output[["deseq2_result_dir"]]
  if (! base::file.exists(output_prefix)) {
    base::dir.create(path = output_prefix,recursive = TRUE)
  }

  # Building command lines for both wald results and fc schinkage
  results_cmd <- base::paste0("DESeq2::results(", extra_results, ")")
  base::message("Command line used for TSV results creation:")
  base::message(results_cmd)

  extra_schrink <- "object = wald, res = results_frame, contrast = contrast"
  if ("extra_schrink" %in% base::names(x = snakemake@params)) {
    extra_schrink <- base::paste(
      ", ",
      base::as.characterx = snakemake@params[["extra_schrink"]]
    )
  }
  schrink_cmd <- base::paste0("DESeq2::lfcSchrink(", extra_schrink, ")")
  base::message("Command line used for log(FC) schrinkage:")
  base::message(schrink_cmd)

  for (resultname in wald_results_names) {
    # Building table
    base::message(base::paste("Saving results for", resultname))
    results_frame <- base::eval(base::parse(text = results_cmd))
    schrink_frame <- base::eval(base::parse(text = schrink_cmd))
    results_frame$log2FoldChange <- schrink_frame$log2FoldChange

    results_path <- base::file.path(
      output_prefix,
      base::paste0("Deseq2_", resultname, ".tsv")
    )

    # Saving table
    utils::write.table(
      x = results_frame,
      file = results_path,
      quote = FALSE,
      sep = "\t",
      row.names = TRUE
    )
  }
}


# If user provides contrasts, then a precise result
# can be extracted from DESeq2 object.
if ("contrast" %in% base::names(snakemake@params)) {
  contrast_length <- base::length(snakemake@params[["contrast"]])
  message(snakemake@params[["contrast"]], contrast_length)

  extra_results <- "object=wald"
  contrast <- NULL

  if (contrast_length == 1) {
    # Case user provided a result name in the `contrast` parameter
    contrast <- base::as.character(x=snakemake@params[["contrast"]])
    contrast <- base::paste0("name='", contrast[1], "'")

  } else if (contrast_length == 2) {
    # Case user provided both tested and reference level
    # In that order! Order matters.
    contrast <- sapply(
      snakemake@params[["contrast"]],
      function(extra) base::as.character(x=extra)
    )
    contrast <- base::paste0(
      "contrast=list('", contrast[1], "', '", contrast[2], "')"
    )

  } else if (contrast_length == 3) {
    # Case user provided both tested and reference level,
    # and studied factor.
    contrast <- sapply(
      snakemake@params[["contrast"]],
      function(extra) base::as.character(x=extra)
    )
    contrast <- base::paste0(
      "contrast=c('",
      contrast[1],
      "', '",
      contrast[2],
      "', '",
      contrast[3],
      "')"
    )

    # Finally saving results as contrast has been
    # built from user input.
    extra_results <- base::paste(extra_results, contrast, sep=", ")
    results_cmd <- base::paste0("DESeq2::results(", extra_results, ")")
    base::message("Result extraction command: ", results_cmd)
    results_frame <- base::eval(base::parse(text = results_cmd))


  }


  if ("filter_theta" %in% base::names(snakemake@output)) {
    base::message("Saving ThetaFilter information: ")
    theta <- metadata(results_frame)$filterNumRej
    # Saving theta filter table on demand
    utils::write.table(
      x = theta,
      file = base::as.character(x = snakemake@output[["filter_theta"]]),
      quote = FALSE,
      sep = "\t",
      row.names = TRUE
    )
  }

  if ("metadata" %in% base::names(snakemake@output)) {
    base::message("Saving Metadata information")
    metadata_table <- data.frame(metadata(results_frame)$filterThreshold)
    metadata_table$filterTheta <- metadata(results_frame)$filterTheta
    metadata_table$alpha <- metadata(results_frame)$alpha
    metadata_table$lfcThreshold <- metadata(results_frame)$lfcThreshold

    # Saving theta filter table on demand
    utils::write.table(
      x = metadata_table,
      file = base::as.character(x = snakemake@output[["metadata"]]),
      quote = FALSE,
      sep = "\t",
      row.names = TRUE
    )
  }

  results_frame$filterThreshold <- results_frame$baseMean > metadata(results_frame)$filterThreshold

  # Saving table
  utils::write.table(
    x = results_frame,
    file = base::as.character(x = snakemake@output[["deseq2_tsv"]]),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE
  )
}



if ("shrinked_wald" %in% base::names(snakemake@output)) {
  contrast_length <- base::length(snakemake@params[["contrast"]])
  message(snakemake@params[["contrast"]], contrast_length)

  extra_results <- 'object=wald, type="apeglm"'
  contrast <- NULL

  if (contrast_length == 1) {
    contrast <- base::as.character(x=snakemake@params[["contrast"]])
    contrast <- base::paste0("coef='", contrast[1], "'")

  } else if (contrast_length == 2) {
    contrast <- sapply(
      snakemake@params[["contrast"]],
      function(extra) base::as.character(x=extra)
    )
    contrast <- base::paste0(
      "contrast=list('", contrast[1], "', '", contrast[2], "')"
    )

  } else if (contrast_length == 3) {
    contrast <- sapply(
      snakemake@params[["contrast"]],
      function(extra) base::as.character(x=extra)
    )
    contrast <- base::paste0(
      "contrast=c('",
      contrast[1],
      "', '",
      contrast[2],
      "', '",
      contrast[3],
      "')"
    )
  }
  extra_results <- base::paste(extra_results, contrast, sep=", ")
  results_cmd <- base::paste0("DESeq2::lfcShrink(", extra_results, ")")
  base::message("Shrunken results extraction command: ", results_cmd)
  results_frame <- base::eval(base::parse(text = results_cmd))
  results_frame$filterThreshold <- results_frame$baseMean > metadata(results_frame)$filterThreshold

  # Saving table
  utils::write.table(
    x = results_frame,
    file = base::as.character(x = snakemake@output[["shrinked_wald"]]),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE
  )
}

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message")
base::sink()