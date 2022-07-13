#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for bcftools annotate"""

from snakemake.shell import shell
from snakemake_wrapper_utils.bcftools import get_bcftools_opts

bcftools_opts = get_bcftools_opts(snakemake)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if "info" in snakemake.input.keys():
    extra += f" --header-lines {snakemake.input.info}"

shell(
    "bcftools annotate {bcftools_opts} "
    "{extra} "
    "--annotations {snakemake.input.annotation} "
    "{snakemake.input.calls} "
    " --output {snakemake.output[0]}"
    "{log}"
)
