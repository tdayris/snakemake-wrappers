#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a box plot of each gene counts
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

from os.path import join
from snakemake.utils import makedirs

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

# Build output directory if necessary
if (outdir := snakemake.output["dir"]) != "":
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

# Stack values in order to plot counts
data = pandas.DataFrame(data.stack())
data.reset_index(inplace=True)
data.columns = ["Target_id", "Sample", "Normalized_Counts"]

data["Condition"] = [
    snakemake.params["condition"][sample]
    for sample in data.Sample
]

logging.debug("Head of post processed dataset:")
logging.debug(data.head())


for target in snakemake.params["targets"]:
    logging.info(f"Working on '{target}'")

    # Define output path
    output_path = join(outdir, f"{target}.png")

    # Subset to searched targets
    tmp = data[data["Target_id"] == target]
    logging.debug("Head of plotted data:")
    logging.debug(tmp.head())

    if len(tmp) == 0:
        msg = f"No target named '{target}' in dataset"
        logging.error(msg)
        raise ValueError(msg)

    # Build plot
    seaborn.set(
        style="ticks",
        palette="pastel"
    )

    seaborn.boxplot(
        x = "Target_id",
        y = "Normalized_Counts",
        hue = "Condition",
        data = tmp
    )

    seaborn.despine(
        offset = 10,
        trim = True
    )

    # Save figure
    matplotlib.pyplot.savefig(
        output_path,
        bbox_inches="tight"
    )
    matplotlib.pyplot.clf()
logging.info("Process over.")
