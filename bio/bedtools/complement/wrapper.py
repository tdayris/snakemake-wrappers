#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for bedtools complement"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

shell(
    "bedtools complement "
    " {extra} "
    " -i {snakemake.input.intervals} "
    " -g {snakemake.input.genome} "
    " > {snakemake.output} "
    " {log} "
)
