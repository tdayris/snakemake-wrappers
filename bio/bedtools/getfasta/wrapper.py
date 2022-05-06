#!/usr/bin/env python3
# coding: utf-8

"""Snakemake wrapper for bedtools getfasta."""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2022, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

out_format = ""
if str(snakemake.output).endswith((".tsv", ".txt")):
    out_format = "-tab"
elif str(snakemake.output).endswith(".bed"):
    out_format = "-bedOut"


shell(
    "bedtools getfasta "
    "-fi {snakemake.input.ref} "
    "-bed {snakemake.input.bed} "
    "-fo {snakemake.output.bed} "
    "{out_format} "
    "{extra} "
    "{log"
)