#!/usr/bin/env python3
# coding: utf-8

"""Annotate annot-sv result with census"""

import pandas
import logging


import os.path as op
import snakemake.utils as smkutils

logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)


def read_census(path: str) -> pandas.DataFrame:
    """Load census file in memory"""
    logging.info("Loading Census data")
    census = pandas.read_csv(path, sep=",", index_col=None, header=0)
    logging.debug(census.head())
    return census[["Tier", "Role in Cancer", "Gene Symbol"]]


def read_annot_sv(path: str) -> pandas.DataFrame:
    """Load annot-sv in memory"""
    logging.info("Loading AnnotSV data")
    annotsv = pandas.read_csv(path, sep="\t", index_col=None, header=0)
    logging.debug(annotsv.head())
    return annotsv[
        [
            "AnnotSV_ID",
            "SV_chrom",
            "SV_start",
            "SV_end",
            "SV_length",
            "SV_type",
            "Samples_ID",
            "ID",
            "FILTER",
            "INFO",
            "AnnotSV_ranking_score",
            "CytoBand",
            "Gene_name",
            "Gene_count",
            "ENCODE_blacklist_left",
        ]
    ]


# Loading datasets
annotsv = read_annot_sv(snakemake.input["annot_sv"])
census = read_census(snakemake.input["census"])
logging.info("Dataset loaded")

# Merging splitted lines
result = pandas.merge(
    left=annotsv,
    right=census,
    how="left",
    left_on="Gene_name",
    right_on="Gene Symbol",
    sort=False,
)
logging.info("Merge succeded")
logging.debug(result.head())

dirpath = op.dirname(snakemake.output["census"])
if not op.exists(dirpath):
    smkutils.makdirs(dirpath)

result.to_csv(snakemake.output["census"], sep="\t", header=True, index=False)
logging.info("Process over")
