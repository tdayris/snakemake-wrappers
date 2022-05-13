#!/usr/bin/env python3
# coding: utf-8

"""Snakemake wrapper for bedtools maskfasta"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2022, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "bedtools maskfasta "
    "-fi {snakemake.input.fasta} "
    "-bed {snakemake.input.bed} "
    "-fo {snakemake.output[0]} "
    "{log}"
)
