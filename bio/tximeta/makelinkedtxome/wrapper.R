# __author__ = "Thibault Dayris"
# __copyright__ = "Copyright 2023, Thibault Dayris"
# __email__ = "thibault.dayris@gustaveroussy.fr"
# __license__ = "MIT"

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file <- base::file(snakemake@log[[1]], open = "wt")
base::sink(log.file)
base::sink(log.file, type = "message")

base::library(package = "BiocFileCache", character.only = TRUE)
base::library(package = "tximeta", character.only = TRUE)

# Overload BiocFileCache directory in order to control
# the produced results. This is done is and only if
# user provides "cache" as an output directory.
if ("cache" %in% base::names(snakemake@output)) {
    origin_path <- BiocFileCache::getBFCOption("CACHE")
    new_cache <- base::as.character(x = snakemake@output[["cache"]])

    if (! base::dir.exists(new_cache)) {
        base::dir.create(new_cache)
    }

    BiocFileCache::setBFCOption(arg = "CACHE", value = new_cache)
    base::message(
        "BiocFileCache directory has been set from `",
        origin_path, "` to `", new_cache, "`."
    )
}

# Salmon index can be either a path to a directory
# or multiple paths to salmon index files in a common
# directory.
salmon_index <- snakemake@input[["salmon_index"]]
if (base::inherits(x = salmon_index, what = "character")) {
    # Then user input points to a single file or directory.
    if (base::file.exists(salmon_index) && ! base::dir.exists(salmon_index)) {
        # Then user input points to a file. We need the index *directory*
        salmon_index <- base::dirname(salmon_index)
    } else {
        # Then user input points to a directory. This is good.
        salmon_index <- base::as.character(x = salmon_index)
    }
} else {
    # Then user input points to a collection of paths.
    # Common prefix should be the salmon index directory.

    # Sort path vector, so the most divergent paths are apart from each other
    # This reduces the number of comparisons in order to find the longest
    # common prefix.
    salmon_index <- base::sort(salmon_index)

    # Split first and last element by character in order to compare
    # them two by two.
    path_split <- base::strsplit(
        salmon_index[c(1, base::length(salmon_index))],
        ""
    )

    # Compare character two by two, return their position
    common_chars <- base::match(
        FALSE,
        base::do.call("==", path_split)
    )

    # -1 is here to account for 1-based vectors in R
    common_chars <- common_chars - 1

    if (common_chars == 0) {
        # Then there is no common substring. Index is in $PWD
        # This is the only possible reason, since Sankemake
        # already tests for file existence
        salmon_index <- "."
    } else {
        # Then there is a common substring and we shall
        # test for its existence in the file system.
        salmon_index <- base::substr(salmon_index[1], 1, common_chars)
        if (! base::dir.exists(salmon_index)) {
            base::stop(
                "Could not find a common prefix for salmon directory. ",
                "Please give the path to the main directory as a string ",
                "or let all the file of the salmon index to be in the ",
                "same parent directory."
            )
        }
    }
}

# Finally build linked txome
tximeta::makeLinkedTxome(
    indexDir = salmon_index,
    source = base::as.character(x = snakemake@params[["source"]]),
    organism = base::as.character(x = snakemake@params[["organism"]]),
    release = base::as.character(x = snakemake@params[["release"]]),
    genome = base::as.character(x = snakemake@params[["genome"]]),
    fasta = snakemake@input[["fasta"]],
    gtf = snakemake@input[["gtf"]],
    write = TRUE,
    jsonFile = snakemake@output[["json"]]
)

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
base::sink(type = "message")
base::sink()