#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for bedtools maskfasta"""

from snakemake.shell import shell

## Extract arguments
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "(bedtools maskfasta"
    " {extra}"
    " -fi {snakemake.input.fasta}"
    " -bed {snakemake.input.bed}"
    " -fo {snakemake.output})"
    " {log}"
)
