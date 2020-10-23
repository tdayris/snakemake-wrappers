"""Snakemake wrapper for bioinfokit Volcanoplot"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020"
__email__ = "Thibault Dayris"
__license__ = "MIT"

from bioinfokit import visuz
import csv
import os
import pandas


# Recover both figure format and name from output file name
figname, figtype = os.path.splitext(snakemake.output[0])
figtype = figtype.strip(".").lower()

# Detect delimiter with python
with open(snakemake.input[0], "r") as table:
    dialect = csv.Sniffer().sniff(table.readline())

# Loading dataset
df = pandas.read_csv(
    snakemake.input[0],
    sep=dialect.delimiter,
    **snakemake.params.get("read_csv", {})
)

# Building and saving volcanoplot
visuz.gene_exp.ma(
    df=df,
    show=False,
    figtype=figtype,
    figname=figname,
    **snakemake.params.get("volcano", {})
)
