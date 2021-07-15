## This script generates a table containing GIS (genomic instability scores) from an EaCoN TCN output using the ASCAT segmenter.

## Requirements :
## - CRAN : optparse, devtools, foreach, dplyr
## - GITHUB : https://github.com/aoumess/chromosomes

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file <- file(
  base::as.character(x=snakemake@log[[1]]),
  open="wt"
);
base::sink(log.file);
base::sink(log.file,type="message");

base::library(package="ASCAT", quietly=TRUE);
base::library(package="EaCoN", quietly=TRUE);
base::library(package="chromosomes", quietly=TRUE);

indir <- base::as.character(x=snakemake@params[["indir"]]);

genome <- "hg19";
if ("genome" %in% base::names(snakemake@params)) {
  genome <- base::as.character(x=snakemake@params[["genome"]]);
}

segmenter <- "ASCAT";
if ("segmenter" %in% base::names(snakemake@params)) {
  segmenter <- base::as.character(x=snakemake@params[["segmenter"]]);
}

## Loading chromosomes data
utils::data(list = genome, package = "chromosomes", envir = environment());


## Loading the gammaEval file
samplename <- base::basename(path=indir);
if ("sample_name" %in% base::names(snakemake@params)) {
  samplename <- base::as.character(x=snakemake@params[["samplename"]]);
}
tcnsegdir <- base::paste0(indir, '/', segmenter, '/ASCN/');
gE <- utils::read.table(
    base::paste0(tcnsegdir, samplename, '.gammaEval.txt'),
    header = TRUE,
    sep = "\t"
);


## Looking for the optimal model
gE$round.psi <- base::round(gE$psi)
gE$round.psi.diff <- base::abs(gE$psi - gE$round.psi)
gE$gof.rank <- base::rank(gE$GoF)
gE$rpd.rank <- base::rank(gE$round.psi.diff)
gE$score <- gE$gof.rank - gE$rpd.rank

best.gamma <- base::sprintf("%.2f", gE$gamma[base::which.max(gE$score)])
print(gE)
print(gE[gE$gamma == best.gamma, ]);

## Loading TCN-CBS file
as.res <- base::readRDS(
  file = base::paste0(
    tcnsegdir, '/gamma', best.gamma, '/', samplename, '.ASCN.ASCAT.RDS'
  )
);

## Computing needed values
as.res$segments$TCN <- as.res$segments$nMajor + as.res$segments$nMajor
as.res$segments$Width <- as.res$segments$endpos - as.res$segments$startpos + 1

## Computing scores
ploidy <- as.res$ploidy$ascat
gof <- as.res$goodnessOfFit
### Sum of segments lengths multiplied by the absolute difference to diploidy,
### divided by the genome length.
SKOR1 <- base::sum(
  base::abs(as.res$segments$TCN - 2
) * as.res$segments$Width) / cs$genome.length
### Sum of segments lengths multiplied by the absolute difference
### to exact ASCAT2-predicted ploidy, divided by the genome length.
SKOR2 <- base::sum(
  base::abs(as.res$segments$TCN - ploidy
) * as.res$segments$Width) / cs$genome.length
### Sum of segments lengths multiplied by the absolute difference to rounded ASCAT2-predicted ploidy, divided by the genome length.
SKOR3 <- base::sum(
  base::abs(as.res$segments$TCN - base::round(ploidy)
) * as.res$segments$Width) / cs$genome.length
### Nb of segments
NBSEG <- base::nrow(as.res$segments)
### Compute a score on a per-chromosome basis, taking the longest ACN (in bp)
### as normal basis
library(foreach)
suppressPackageStartupMessages(library(dplyr))
LOKAL1 <- base::sum(foreach (k = unique(as.res$segments$chr), .combine = "c") %do% {
  miniseg <- tibble::as_tibble(as.res$segments[as.res$segments$chr == k,]);
  basistbl <- miniseg %>% dplyr::group_by(
    TCN
  ) %>% dplyr::summarise(base::sum(Width));
  basis <- basistbl$TCN[base::which.max(as.data.frame(basistbl)[,2])];
  KSCOR <- base::sum(base::abs(miniseg$TCN  - basis) * miniseg$Width) / cs$genome.length
  return(KSCOR)
})

## Generating output
best.gamma = base::as.numeric(best.gamma);
outdf <- data.frame(
  best.gamma = best.gamma,
  psi = gE[gE$gamma == best.gamma, "psi"],
  gof = gE[gE$gamma == best.gamma, "GoF"],
  score = gE[gE$gamma == best.gamma, "score"],
  SCORE1 = SKOR1,
  SCORE2 = SKOR2,
  SCORE3 = SKOR3,
  NBSEG = NBSEG,
  LOCAL1 = LOKAL1,
  stringsAsFactors = FALSE
);
utils::write.table(
    x = outdf,
    file = base::paste0(
      indir, '/', samplename, '_GIS_from_best_gamma.txt'
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
);

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type="message");
base::sink();
