#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for WisecondorX convert"""

from snakemake.shell import shell

log = snakemake.log_fmt_shell()
extra = snakemake.params.get("extra", "")

shell(
    "WisecondorX convert {snakemake.input.aln} "
    "--reference {snakemake.input.ref} "
    "{snakemake.output[0]} {extra} {log}"
)
