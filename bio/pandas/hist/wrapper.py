#!/usr/bin/python3.8
# conding: utf-8

"""
Plot a histogram based on a tsv table
"""

import pandas
import matplotlib.pyplot as plt

data = pandas.read_csv(
    snakemake.input["tsv"],
    sep="\t",
    header=0,
    index_col=None
)

data.hist(**snakemake.params["extra"])
plt.savefig(snakemake.output["png"])
