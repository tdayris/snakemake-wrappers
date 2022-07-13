#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for wisecondorx newref"""

from snakemake.shell import shell

log = snakemake.log_fmt_shell()
extra = snakemake.params.get("extra", "")

shell(
    "WisecondorX newref {snakemake.input.aln} "
    "{snakemake.output[0]} {extra} {log}"
)
