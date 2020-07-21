#!/usr/bin/R

# Plot mutational spectrum
library(package = "SomaticSignatures", quietly=TRUE);


contexts <- base::readRDS(
  file = base::as.character(
    x = snakemake@input[['context']]
  )
);

extra <- "contexts";
if ("extra" %in% names(snakemake@params)) {
  extra <- base::paste(
    extra,
    base::as.character(x = snakemake@params[["extra"]]),
    sep = ", "
  );
}

command <- paste0(
  "SomaticSignatures::plotMutationSpectrum(",
  extra,
  ")"
);
print(command);


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
