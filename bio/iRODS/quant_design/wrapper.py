#!/usr/bin/python3.8
# conding: utf-8

"""
Build a rna-count-salmon design file from a TSV formatted iRODS metadata
"""

from snakemake.shell import shell

import pandas

#  Load input table
metadata = pandas.read_csv(
    snakemake.input[0],
    sep="\t",
    header=0,
    index_col=None,
    dtype=str
)


# Merge upstream and downstream files if paired dataset
if "paired" not in snakemake.params.keys() or snakemake.params["paired"]:
    common_keys = snakemake.params["common"]
    aggregate = {
        k: ",".join
        for k in set(metadata.columns) - set(common_keys)
    }

    metadata.groupby(common_keys).agg(aggregate)

#  Optional column rename to fit rna-count-salmon requirements
if "rename" in snakemake.params.keys():
    metadata.rename(snakemake.params["rename"])


# Save output files
metadata.to_csv(
    snakemake.output[0],
    sep="\t",
    header=True,
    index=True
)

shell("sed -i 's/,/\\t/g' {snakemake.output[0]}")
