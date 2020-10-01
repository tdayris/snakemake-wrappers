#!/usr/bin/R

# This script downloads and build mouse to human correspondancy

# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2020, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Loading libraries
base::library(package = "biomaRt", quietly = TRUE);
base::library(package = "dplyr", quietly = TRUE);
base::message("Library loaded");

# Loading gene list
mouse_table <- utils::read.table(
  file=snakemake@input[["mouse"]],
  header=TRUE,
  sep="\t",
  stringsAsFactors=FALSE
);

mgi_id <- "gene_id";
if ("mgi_id" %in% base::names(snakemake@params)) {
  mgi_id <- base::as.character(
    x=snakemake@params[["mgi_id"]]
  );
}

mouse_genes <- mouse_table[, mgi_id];

# Loading annotations from ensembl
human <- biomaRt::useMart(
  "ensembl",
  dataset="hsapiens_gene_ensembl",
);

mouse <- biomaRt::useMart(
  "ensembl",
  dataset="mmusculus_gene_ensembl"
);
base::message("Marts retrieved");

genes <- biomaRt::getLDS(
  attributes = c("mgi_symbol"),
  filters = "mgi_symbol",
  values = mouse_genes,
  mart = mouse,
  attributesL = c("hgnc_symbol"),
  martL = human,
  uniqueRows = TRUE
);
base::message("Namespace merged");

mousexhuman <- unique(genes[, 2]);
nb_mouse_genes <- length(mouse_genes);
nb_human_genes <- length(mousexhuman);

if (nb_mouse_genes != nb_human_genes) {
  base::message(
    base::paste(
      "Some genes were lost:",
      nb_human_genes,
      "/",
      nb_mouse_genes,
      "=",
      (nb_human_genes/nb_mouse_genes)*100,
      "%"
    )
  );
  base::print(setdiff(mouse_genes, genes$HGNC.symbol))

} else {
  base::message("All genes were translated");
}

results <- base::merge(
  x = mouse_table,
  y = genes,
  by.x = mgi_id,
  by.y = "MGI.symbol",
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
);

results[, mgi_id] <- NULL;

results %>% dplyr::group_by(
  HGNC.symbol
) %>% dplyr::summarise_if(
  base::is.numeric,
  base::sum,
  na.rm=TRUE
) %>% dplyr::filter(
  !is.na(HGNC.symbol)
) -> results;

utils::write.table(
  x=results,
  file=snakemake@output[["translated"]],
  sep="\t",
  quote=FALSE,
  col.names=TRUE,
  row.names=FALSE
);
