#!/usr/bin/python3.8
# conding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

"""
Plot the distribution
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
    snakemake.input["tsv"],
    sep="\t",
    header=0
)
logging.debug("Dataframe:")
logging.debug(data.head())


seaborn.displot(
    data=data,
    **snakemake.params.get("extra", {})
)
logging.debug("Plot created")

matplotlib.pyplot.savefig(
    snakemake.output["png"],
    bbox_inches="tight"
)
logging.debug("Process over.")
