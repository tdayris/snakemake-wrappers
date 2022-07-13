#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2021, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# This script takes a count matrix, a coldata and a formula and returns
# a rds-formatted dds object.

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt");
base::sink(log.file);
base::sink(log.file,type="message");

base::library(package="DESeq2", quietly=TRUE);

counts <- utils::read.table(
  file=snakemake@input[["counts"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
);

gene_names_col <- "genes";
if ("count_gene_col" %in% base::names(snakemake@params)) {
  gene_names_col <- base::as.character(x=snakemake@params[["count_gene_col"]]);
}
rownames(counts) <- counts[, gene_names_col];

numerical_columns <- base::unlist(base::lapply(x, base::is.numeric));
counts <- counts[, numerical_columns];

count_filter <- 0;
if ("count_filter" %in% names(snakemake@params)) {
  count_filter <- base::as.numeric(
    x = snakemake@params[["count_filter"]]
  );
}

keep <- rowSums(counts(counts)) > count_filter;
counts <- counts[keep, ];
base::message("Counts loaded");

coldata <- utils::read.table(
  file=snakemake@input[["coldata"]],
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
);
base::message("Coldata loaded");

extra <- "countData=counts, colData=coldata"
if ("extra" %in% base::names(snakemake@params)) {
  extra <- base::paste(
    extra,
    base::as.character(x=snakemake@params[["extra"]])
    sep=","
  );
}

# Create object
dds <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::DESeqDataSetFromMatrix(", extra, ");"
    )
  )
);

levels <- ("levels" %in% base::names(snakemake@params));
factor <- ("factor" %in% base::names(snakemake@params));
if (levels & factor) {
  # Levels should come with reference as first item
  levels <- sapply(
    snakemake@params[["levels"]],
    function(level) base::as.character(x = level)
  );
  factor <- base::as.character(x = snakemake@params[["factor"]]);

  dds[[factor]] <- base::factor(dds[[factor]], levels = levels);
  dds[[factor]] <- stats::relevel(dds[[factor]], ref = levels[[1]]);
  dds[[factor]] <- droplevels(dds[[factor]]);
  base::write("Factors have been releveled", stderr());
}


# Save as RDS
output_path <- base::as.character(x = snakemake@output[["dds"]]);
base::saveRDS(
  obj = dds,
  file = output_path
);
base::message("Process over");

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
