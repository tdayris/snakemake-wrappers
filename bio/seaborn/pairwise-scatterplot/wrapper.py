#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a clustered heatmap of multiple sample
"""

import logging
import matplotlib
import matplotlib.pyplot
import numpy
import pandas
import seaborn

from os.path import basename, dirname
from snakemake.utils import makedirs

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

# Build output directory if necessary
if (outdir := basename(dirname(snakemake.output["png"]))) != "":
    makedirs(outdir)
    logging.debug(f"Directory: '{outdir}' created.")

data = pandas.read_csv(
    snakemake.input["counts"],
    sep="\t",
    index_col=0,
    header=0
)

logging.debug("Loaded dataset:")
logging.debug(data.head())

seaborn.set(
    style="ticks",
    color_codes=True
)

g = seaborn.pairplot(
    data,
    kind="scatter",
    dropna=True,
    diag_kind="kde"
)

matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)
logging.info("Process over")
