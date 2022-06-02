#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for WisecondorX predict"""

from snakemake.shell import shell

log = snakemake.log_fmt_shell()
extra = snakemake.params.get("extra", "")
prefix = snakemake.params.get("prefix", None)
if prefix is None:
    prefix = snakemake.output["bins"][:-len("_bins.bed")]

shell(
    "WisecondorX predict {snakemake.input.test} "
    "{snakemake.input.ref} {prefix} {extra} {log}"
)
