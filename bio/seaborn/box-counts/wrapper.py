#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a box plot of each sample counts
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

# Load normalized counts
data = pandas.read_csv(
    snakemake.input["counts"],
    sep="\t",
    header=0,
    index_col=0
)

# Remove possible text annotations and validate
data = data[list(data.select_dtypes(include=[numpy.number]).columns.values)]

drop_null = snakemake.params.get("drop_null", False)
if drop_null is True:
    data = data.loc[~(data == 0).all(axis=1)]
    data.dropna(axis=0, how="all", inplace=True)

# Stack values in order to plot counts
data = pandas.DataFrame(data.stack())
data.reset_index(inplace=True)
data.columns = ["Target_id", "Sample", "Normalized_Counts"]

logging.debug("Head of the post-processed dataframe:")
logging.debug(data.head())

seaborn.set(
    style="ticks",
    palette="pastel"
)

seaborn.boxplot(
    y = "Sample",
    x = "Normalized_Counts",
    data = data,
    orient="h"
)

seaborn.despine(
    offset = 10,
    trim = True
)

matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)
logging.debug("Process over.")
