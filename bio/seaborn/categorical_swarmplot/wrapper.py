#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a categorical swarmplot over dataframes
"""


__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import matplotlib
import matplotlib.pyplot
import pandas
import seaborn

import os.path

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

seaborn.set(style="white", context="talk")

categories = pandas.DataFrame.from_dict(
    data={
        os.path.splitext(os.path.basename(frame))[0]: pandas.read_csv(
            frame,
            sep=("\t" if frame.endswith(".tsv") else ","),
            header=0,
            index_col=None
        )[snakemake.params["column"]].value_counts(sort=True).to_dict()
        for frame in snakemake.input
    },
    orient="index"
)

logging.debug("Frequency table:")
logging.debug(categories.head())

if "tsv" in snakemake.output.keys():
    categories.to_csv(
        snakemake.output["tsv"],
        sep="\t",
        header=True,
        index=True
    )

if "png" in snakemake.output.keys():
    categories = categories.reset_index().melt(
        id_vars="index",
        var_name="Categories",
        value_name="Frequency"
    )

    seaborn.catplot(
        x = "Categories",
        y = "Frequency",
        hue = "index",
        data = categories
    )

    matplotlib.pyplot.savefig(
        snakemake.output["png"],
        bbox_inches="tight"
    )
