#!/usr/bin/python3.8
# -*- coding: utf-8 -*-

"""
split complex differential analysis design into simple ones
"""

import os
import csv
import logging
import pandas

from itertools import combinations
from snakemake.utils import makedirs


# Initiate logging behaviour
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)


# Build output directory
out_design_dir = snakemake.output.get("out_dir", None)
if out_design_dir is None:
    out_design_dir = os.path.dirname(snakemake.output[0])
makedirs(out_design_dir)
logging.info("Output directory built")


# Loading input data
with open(snakemake.input["design"], "r") as table:
    dialect = csv.Sniffer().sniff(table.readline())
    logging.debug(f"Detected dialect: {dialect.delimiter}")

read_csv_params = dict(
    sep=dialect.delimiter,
    header=0,
    index_col=0
)
if "read_csv" in snakemake.params.keys():
    read_csv_params.update(**snakemake.params["read_csv"])

complete_design = pandas.read_csv(
    snakemake.input["design"],
    **read_csv_params
)
logging.info("Input design loaded")

if "columns_to_aggregate" in snakemake.params.keys() and snakemake.params["columns_to_aggregate"] is not None:
    for cols in snakemake.params["columns_to_aggregate"]:
        complete_design["_".join(cols)] = (
            complete_design[cols].astype(str)
                                 .apply("_".join, axis=1)
                                 .str.replace(" ", "_")
                                 .str.strip()
        )
        logging.info("%s aggregated", str(cols))

if "columns_to_remove" in snakemake.params.keys() and snakemake.params["columns_to_remove"] is not None:
    removed = snakemake.params["columns_to_remove"]
    complete_design.drop(removed, axis=1, inplace=True)
    logging.info("Columns %s removed", str(removed))

for col in complete_design.columns:
    logging.info("Considering factor %s", col)
    levels = sorted(
        k for k, v in complete_design[col].value_counts().items() if v > 1
    )

    if len(levels) <= 1:
        # We need at least two levels to perfom a differential analysis
        continue


    for l1, l2 in combinations(levels, 2):
        # R and DESeq2 sort levels through alphanumerical order. We have
        # to provide user-understandable name ; so let us sort out the
        # levels and guess reference name.
        for level, ref in [[l1, l2], [l2, l1]]:

            # Building humand readable design name
            design_name = f"{out_design_dir}/DGE_considering_factor_{col}_comparing_test_{level}_vs_reference_{ref}.tsv"
            tmp = complete_design[complete_design[col].isin([level, ref])]
            tmp.replace([" ", level, ref], ["", f"test_{level}", f"reference_{ref}"], inplace=True)
    
            logging.info(
                "%s: Considering level %s and reference %s",
                design_name,
                level,
                ref
            )
    
            tmp.to_csv(
                design_name,
                sep="\t",
                index=True,
                header=True,
                index_label="Sample_id"
            )