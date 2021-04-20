#!/usr/bin/python

"""
Snakemake wrapper for GATK CollectSequencingArtifactMetrics
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

dbsnp = ""
if "dbsnp" in snakemake.input.keys():
    dbsnp = "--DB_SNP {}".format(snakemake.input["dbsnp"])

intervals = ""
if "intervals" in snakemake.input.keys():
    intervals = "--INTERVALS {}".format(snakemake.input["intervals"])

prefix = snakemake.params.get("prefix", "")
if prefix == "":
    prefix = snakemake.output[0][:-len(".pre_adapter_detail_metrics")]


shell(
    "gatk --java-options '{java_opts}' "  # Base GATK Call with java arguments
    "CollectSequencingArtifactMetrics {extra} "  # Optional arguments
    "{dbsnp} "  # Optional dbsnp file
    "{intervals} "  # Optional bed file
    "--INPUT {snakemake.input.bam} "  # Path to input bam file
    "--REFERENCE_SEQUENCE {snakemake.input.ref} "  # Path to reference genome
    "--OUTPUT {prefix} {log}"  # Path to output file
)
