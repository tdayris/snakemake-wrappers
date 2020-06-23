#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for deeptools computeMatrix"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

subcommand = snakemake.params.get("subcommand", "scale-regions")

extra = snakemake.params.get("extra", "")
if "blacklist" in snakemake.input.keys():
    extra += "--blackListFileName {snakemake.input.blacklist}"

output = ""
if "gzip" in snakemake.output.keys():
    output += " --outFileName {} ".format(snakemake.output.gzip)
if "matrix" in snakemake.output.keys():
    output += " --outFileNameMatrix {} ".format(snakemake.output.matrix)
if "regions" in snakemake.output.keys():
    output += " --outFileSortedRegions {} ".format(snakemake.output.regions)
if len(output) == 0:
    raise KeyError(
        "Output section should include: gzip, matrix, and/or regions"
    )

shell(
    "computeMatrix "
    " {subcommand} "
    " --regionsFileName {snakemake.input.regions} "
    " --scoreFileName {snakemake.input.scores} "
    " --numberOfProcessors {snakemake.threads} "
    " {output} "
    " {extra} "
    " {log} "
)
