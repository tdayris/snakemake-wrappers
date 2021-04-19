#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for bcftools query"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")
output = snakemake.output[0]
call = snakemake.input[0]

if "sample_list" in snakemake.output.keys():
    extra = "-l"
    output = snakemake.output["sample_list"]

elif "bed" in snakemake.output.keys():
    extra = "-f'%CHROM\t%POS0\t%END\t%ID\n'"
    output = snakemake.output["bed"]

shell(
    "bcftools query "
    "{extra} "
    "{call} "
    "> {output} "
    "{log}"
)
