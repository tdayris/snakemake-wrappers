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

intervals = snakemake.input.intervals
if isinstance(intervals, str):
    intervals = [intervals]
intervals = list(map("--intervals {}".format, intervals))

bams = snakemake.input.bam
if isinstance(bams, str):
    bams = [bams]
bams = list(map("--input {}".format, bams))

variants = "--variant {}".format(snakemake.input["vcf"])

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)


shell(
    "gatk --java-options '{java_opts}' GetPileupSummaries {extra} "
    "{bams} --output {snakemake.output.table} {log}"
)
