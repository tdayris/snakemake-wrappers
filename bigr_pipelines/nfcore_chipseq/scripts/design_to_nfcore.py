#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script creates a design.csv usable by nextflow
"""

import pandas
import logging

logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    level=logging.DEBUG
)

design = snakemake.params["design"]
fastq_links = snakemake.params["fastq_links"]

nxf_list = [
    {"group": sample, "fastq_1": fastq["r1"], "fastq_2": fastq["r2"]} 
    for sample, fastq in fastq_links.items()
]
nxf_design = pandas.DataFrame.from_records(nxf_list)

del design["Upstream_file"]
del design["Downstream_file"]

nxf_design = pandas.merge(
    left=nxf_design,
    right=design,
    left_on="group",
    right_on="Sample_id",
    how="left"
)

nxf_design.to_csv(
    snakemake.output[0],
    sep=",",
    index=False
)