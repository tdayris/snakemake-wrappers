"""Snakemake wrapper for bioinfokit HeatMap"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020"
__email__ = "Thibault Dayris"
__license__ = "MIT"

from bioinfokit import visuz
import csv
import logging
import os
import pandas


logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)
logging.getLogger('matplotlib.font_manager').disabled = True


# Recover both figure format and name from output file name
figname, figtype = os.path.splitext(snakemake.output[0])
figtype = figtype.strip(".").lower()

# Detect delimiter with python
with open(snakemake.input[0], "r") as table:
    dialect = csv.Sniffer().sniff(table.readline())
    logging.debug(f"Detected dialect: {dialect.delimiter}")

# Loading dataset
df = pandas.read_csv(
    snakemake.input[0],
    sep=dialect.delimiter,
    **snakemake.params.get("read_csv", {})
)
logging.debug(df.head())

# Building and saving volcanoplot
visuz.gene_exp.hmap(
    df=df,
    show=False,
    figtype=figtype,
    figname=figname,
    **snakemake.params.get("hmap", {})
)
