#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for bedtools getfasta"""

from snakemake.shell import shell

## Extract arguments
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if str(snakemake.output[0]).endswith("bed"):
    extra += " -bedOut"
elif str(snakemake.output[0]).endswith("tsv"):
    extra += " -tab"


shell(
    "(bedtools getfasta"
    " {extra}"
    " -fi {snakemake.input.fasta}"
    " -bed {snakemake.input.bed}"
    " -fo {snakemake.output})"
    " {log}"
)
