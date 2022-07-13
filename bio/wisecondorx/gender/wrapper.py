#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for WisecondorX gender"""

from snakemake.shell import shell

log = snakemake.log_fmt_shell()
extra = snakemake.params.get("extra", "")

shell(
    "WisecondorX gender {snakemake.input.test} "
    "{snakemake.input.ref} {extra} > {snakemake.output} "
    "{log}"
)
