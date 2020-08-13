#!/usr/bin/R

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

library(package = "pheatmap", quietly = TRUE);
library(package = "DESeq2", quietly = TRUE);

# Loading datasets and parameters
dds <- base::readRDS(
  file = base::as.character(x = snakemake@input[["dds"]])
);

norm <- base::readRDS(
  file = base::as.character(x = snakemake@input[["norm"]])
);

condition_array <- sapply(
  snakemake@params[["condition_array"]],
  function(cond) base::as.character(x = cond)
);

padj_threshold <- 0.05;
if ("padj_threshold" %in% names(snakemake@params)) {
  padj_threshold <- snakemake@params[["padj_threshold"]];
}

annotation <- base::as.data.frame(
  x = colData(dds)[, condition_array]
);

# Buildcommand line
extra <- 'normalized_counts, annotation_col = annotation';
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
    extra,
    snakemake@params[["extra"]],
    sep=", "
  );
}

command <- paste0(
  "pheatmap::pheatmap(",
  extra,
  ")"
);
print(command);


# For each result in this deseq2 object, print pheatmap
names <- DESeq2::resultsNames(object = dds);

for (resultname in names) {
  base::message(base::paste("Building heatmap for", resultname));
  results_frame <- DESeq2::results(
    object = dds,
    name = resultname,
    independentFiltering = TRUE,
    # alpha = padj_threshold,
    pAdjustMethod = "BH",
    cooksCutoff = TRUE
  );

  not_na <- complete.cases(results_frame);
  significant <- results_frame$padj <= padj_threshold;
  select_ids <- row.names(results_frame[not_na & significant, ]);

  normalized_counts <- assay(norm)[select_ids, ];
  print(head(assay(norm)))
  print(dim(assay(norm)))
  print(head(select_ids))
  print(head(normalized_counts))
  print(dim(normalized_counts))

  # Build plot and save it
  png(
    filename = snakemake@output[["png"]],
    width = 1024,
    height = 768,
    units = "px",
    type = "cairo"
  );

  base::eval(
    base::parse(
      text = command
    )
  );

  dev.off();
}
