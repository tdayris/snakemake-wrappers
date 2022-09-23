#!/usr/bin/python3
# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2022, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake_wrapper_utils.bcftools import get_bcftools_opts

bcftools_opts = get_bcftools_opts(snakemake)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "bcftools stats {bcftools_opts} "
    "{extra} "
    "{snakemake.input[0]} "
    "-o {snakemake.output[0]} "
    "{log}"
)