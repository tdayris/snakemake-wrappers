#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for deeptools bamcoverage"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")
if "regions" in snakemake.wildcards.keys():
    extra += " --region {} ".format(snakemake.wildcards["regions"])


if "blacklist" in snakemake.input.keys():
    extra += " --blackListFileName {} ".format(snakemake.input["blacklist"])

output = snakemake.output["coverage"]
if output.endswith(".bw"):
    extra += " --outFileFormat bigwig "
else:
    extra += " --outFileFormat bedgraph "

shell(
    "bamCoverage "
    " --numberOfProcessors {snakemake.threads} "
    " {extra} "
    " --bam {snakemake.input.bam} "
    " --outFileName {output}"
    " {log} "
)
