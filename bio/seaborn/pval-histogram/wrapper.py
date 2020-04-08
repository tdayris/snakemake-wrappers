#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a pvalue histogram
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn

from os.path import join, basename, dirname
from snakemake.utils import makedirs

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

# Build output directory if necessary
if (outdir := basename(dirname(snakemake.output["png"]))) != "":
    makedirs(outdir)
    logging.debug(f"Directory '{outdir}' created")


# Load dataset
data = pandas.read_csv(
    snakemake.input["deseq2"],
    sep="\t",
    header=0,
    index_col=0
)
data = data["padj"]

# Build intervals
data = data.value_counts(
    bins = 20,
    dropna = True,
    ascending = False,
    sort = False
)
logging.debug("Head of pval bins:")
logging.debug(data.head())

# Build graph
seaborn.set(
    style="darkgrid",
    color_codes="muted"
)
f = seaborn.barplot(
    x = data.index,
    y = data.values,
    color = "b"
)

f.axes.get_xaxis().set_visible(False)


matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)
logging.info("Process over")
