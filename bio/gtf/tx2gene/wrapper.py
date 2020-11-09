#!/usr/bin/python3.8
# conding: utf-8

"""
This script extracts a tsv from a GTF. This tsv contains gene_id,
transcript_id, gene_name. Optionaly, headers and positions can
also be given.
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import numpy
import pandas


logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


# Read GTF file
gtf = pandas.read_csv(
    snakemake.input["gtf"],
    sep="\t",
    header = None,
    comment = "#"
)
gtf.columns = [
    "chromosome", "source", "feature", "start", "stop",
    "confidence", "strand", "frame", "annotation"
]
gtf.set_index(
    ["chromosome", "source", "strand", "feature"],
    inplace = True
)
logging.debug("Head of the GTF file:")
logging.debug(gtf.head())
logging.info("GTF loaded")


# Extract transcripts identifiers, gene identifiers and gene names
gtf["gene_id"] = gtf["annotation"].str.extract(r"gene_id \"(\w+)\"")
gtf["transcript_id"] = gtf["annotation"].str.extract(r"transcript_id \"(\w+)\"")
gtf["gene_name"] = gtf["annotation"].str.extract(r"gene_name \"(\w+)\"")
del gtf["annotation"]


# Filter on transcripts to reduce the amound of data
gtf[gtf.index.get_level_values("feature") == "transcript"]
gtf.index = gtf.index.droplevel("feature")

logging.debug("Head of the GTF file:")
logging.debug(gtf.head())
logging.info("GTF parsed")


# Save results and their subsets on demand
if "tx2gene" in snakemake.output.keys():
    tmp = gtf.copy()
    tmp = tmp[["gene_id", "transcript_id", "gene_name"]]
    tmp.to_csv(
        snakemake.output["tx2gene"],
        sep = "\t",
        index = False,
        header = False
    )
    logging.info("Tx2Gene table saved")

if "tx2gene_large" in snakemake.output.keys():
    tmp = gtf.copy()
    tmp = tmp[["gene_id", "transcript_id", "gene_name", "start", "stop"]]
    tmp.to_csv(
        snakemake.output["tx2gene_large"],
        sep = "\t",
        index = True,
        header = True
    )
    logging.info("Tx2Gene table saved with positions and source info")

if "gene2gene" in snakemake.output.keys():
    tmp = gtf.copy()
    tmp = tmp[["gene_id", "gene_name"]]
    tmp.drop_duplicates(inplace=True)
    tmp.to_csv(
        snakemake.output["gene2gene"],
        sep = "\t",
        index = False,
        header = False
    )
    logging.info("Gene2Gene table saved")


if "gene2gene_large" in snakemake.output.keys():
    tmp = gtf.copy()
    tmp = tmp[["gene_id", "gene_name"]]
    tmp.drop_duplicates(inplace=True)
    tmp.to_csv(
        snakemake.output["gene2gene_large"],
        sep = "\t",
        index = True,
        header = True
    )
    logging.info("Gene2Gene table saved with chromosome and source info")
