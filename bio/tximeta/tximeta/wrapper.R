# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]], open = "wt")
base::sink(log.file)
base::sink(log.file, type = "message")


base::library(package = "tximeta", character.only = TRUE)


tximeta::loadLinkedTxome(jsonFile = snakemake@input[["jsonfile"]])

coldata <- utils::read.table(
    file = snakemake@input[["coldata"]],
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)
if (! "files" %in% base::colnames(coldata)) {
    base::stop(
        "`file` column should be in coldata. This column ",
        "describes the path to the salmon quant file for ",
        "a given sample."
    )
} else if (! "names" %in% base::colnames(coldata)) {
   base::stop(
    "`names` column should be in coldata. This column ",
    "describes the name of a given sample."
   )
} else {
    base::message("The coldata contains required columns.")
}

extra_tximeta <- "coldata = coldata"
if ("extra_tximeta" %in% base::names(snakemake@params)) {
    extra_tximeta <- base::paste(
        extra_tximeta,
        base::as.character(x = snakemake@params[["extra_tximeta"]]),
        sep = ", "
    )
}

tximeta_cmd <- base::paste0(
    "tximeta::tximeta(",
    extra_tximeta,
    ")"
)

txi <- base::eval(base::parse(text = command))

base::saveRDS(
    object = txi,
    file = snakemake@output[["txi"]]
)

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()