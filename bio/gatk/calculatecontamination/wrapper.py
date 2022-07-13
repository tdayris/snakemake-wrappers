#!/usr/bin/python

"""
Snakemake wrapper for GATK GetPileupSummaries
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import os

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

matched = ""
if "normal" in snakemake.input.keys():
    matched = "--matched-normal {}".format(snakemake.input["normal"])


segmentation = ""
if "segmentation" in snakemake.output.keys():
    segmentation = "--tumor-segmentation {}".format(
        snakemake.output["segmentation"]
    )


shell(
    "gatk --java-options '{java_opts}' CalculateContamination {extra} "
    "--input {snakemake.input.summary} --output {snakemake.output.table} "
    "{matched} {segmentation} {log}"
)
