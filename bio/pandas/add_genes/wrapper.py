#!/usr/bin/python3.8
# conding: utf-8

"""
Join multiple two tsv files on a given column. One of
these tsv file if formatted as follows:

1: gene_id
2: transcript_id
3: gene_name

Optional columns can be:

4: chr
5: start
6: end
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import logging
import pandas

from os.path import basename, dirname
from snakemake.utils import makedirs

if (outdir := dirname(snakemake.output["tsv"])) != "":
    makedirs(outdir)

genes = snakemake.params.get("genes", False)

annot = pandas.read_csv(
    snakemake.input["tx2gene"],
    sep="\t",
    header=(
        0
        if snakemake.params.get("header", None) is not None
        else None
    ),
    index_col=None,
    dtype=str
)
annot.columns = (
    ["Ensembl_ID", "Hugo_ID"]
    if genes is True
    else ["Ensembl_ID", "Transcript_ID", "Hugo_ID"]
)

tsv = pandas.read_csv(
    snakemake.input["tsv"],
    sep="\t",
    header=0,
    index_col=0
)

merged_frame = pandas.merge(
    tsv,
    annot,
    left_index=True,
    right_on=(
        "Ensembl_ID"
        if genes is True
        else "Transcript_ID"
    ),
    how="left"
)

merged_frame.to_csv(
    snakemake.output["tsv"],
    sep="\t",
    index=False,
    header=True
)
