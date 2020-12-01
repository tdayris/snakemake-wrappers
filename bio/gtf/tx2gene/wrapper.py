#!/usr/bin/python3.8
# conding: utf-8

"""
This script extracts a tsv from a GTF. This tsv contains Gene_ID,
transcript_id, Gene_Name. Optionaly, headers and positions can
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


def fill_empty(frame: pandas.DataFrame,
               fill_col: str,
               ref_col: str) -> pandas.DataFrame:
    """
    Fill NA and empty values in a frame
    """
    frame[fill_col] = [
        ref if ((fill is None) or (fill == numpy.NaN)) else fill
        for fill, ref in zip(frame[fill_col], frame[ref_col])
    ]
    return frame


# Read GTF file
gtf = pandas.read_csv(
    snakemake.input["gtf"],
    sep="\t",
    header = None,
    comment = "#"
)
gtf.columns = [
    "Chromosome", "source", "feature", "Start", "Stop",
    "confidence", "Strand", "frame", "annotation"
]
gtf.set_index(
    ["Chromosome", "source", "Strand", "feature"],
    inplace = True
)
logging.debug("Head of the GTF file:")
logging.debug(gtf.head())
logging.info("GTF loaded")


# Extract transcripts identifiers, gene identifiers and gene names
gtf["Gene_ID"] = gtf["annotation"].str.extract(r"gene_id \"(\w+)\"")
gtf["transcript_id"] = gtf["annotation"].str.extract(r"transcript_id \"(\w+)\"")
gtf["Gene_Name"] = gtf["annotation"].str.extract(r"gene_name \"(\w+)\"")
del gtf["annotation"]


# Filter on transcripts to reduce the amound of data
gtf[gtf.index.get_level_values("feature") == "transcript"]
gtf.index = gtf.index.droplevel("feature")

logging.debug("Head of the GTF file:")
logging.debug(gtf.head())
logging.info("GTF parsed")

# Removing patches on identifiers
drop_patches = snakemake.params.get("drop_patches", True)
if drop_patches is True:
    gtf["Gene_ID"] = gtf["Gene_ID"].str.split(r"\.", expand=True)[0]
    gtf["transcript_id"] = gtf["transcript_id"].str.split(r"\.", expand=True)[0]


# Save results and their subsets on demand
if "tx2gene_small" in snakemake.output.keys():
    tmp = gtf.copy()
    tmp = tmp[["transcript_id", "Gene_ID"]]
    tmp = fill_empty(tmp.copy(), "Gene_ID", "transcript_id")
    tmp.drop_duplicates(inplace=True)
    tmp.to_csv(
        snakemake.output["tx2gene_small"],
        sep = "\t",
        index = False,
        header = False
    )
    logging.info("Tx2Gene table saved without gene names")

if "tx2gene" in snakemake.output.keys():
    tmp = gtf.copy()
    tmp = tmp[["Gene_ID", "transcript_id", "Gene_Name"]]
    tmp.drop_duplicates(inplace=True)
    tmp = fill_empty(tmp.copy(), "Gene_ID", "transcript_id")
    tmp = fill_empty(tmp.copy(), "Gene_Name", "Gene_ID")
    tmp.to_csv(
        snakemake.output["tx2gene"],
        sep = "\t",
        index = False,
        header = False
    )
    logging.info("Tx2Gene table saved")

if "tx2gene_large" in snakemake.output.keys():
    tmp = gtf.copy()
    tmp = tmp[["Gene_ID", "transcript_id", "Gene_Name", "Start", "Stop"]]
    tmp = fill_empty(tmp.copy(), "Gene_ID", "transcript_id")
    tmp = fill_empty(tmp.copy(), "Gene_Name", "Gene_ID")
    tmp.drop_duplicates(inplace=True)
    tmp.to_csv(
        snakemake.output["tx2gene_large"],
        sep = "\t",
        index = True,
        header = True
    )
    logging.info("Tx2Gene table saved with positions and source info")

if "gene2gene" in snakemake.output.keys():
    tmp = gtf.copy()
    tmp = tmp[["Gene_ID", "Gene_Name"]]
    tmp = fill_empty(tmp.copy(), "Gene_Name", "Gene_ID")
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
    tmp = tmp[["Gene_ID", "Gene_Name"]]
    tmp = fill_empty(tmp.copy(), "Gene_Name", "Gene_ID")
    tmp.drop_duplicates(inplace=True)
    tmp.to_csv(
        snakemake.output["gene2gene_large"],
        sep = "\t",
        index = True,
        header = True
    )
    logging.info("Gene2Gene table saved with Chromosome and source info")
