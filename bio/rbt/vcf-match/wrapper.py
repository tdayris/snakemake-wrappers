#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for rust-bio-tools vcf-match"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

shell(
    "rbt vcf-match {snakemake.input.reference} {extra} "
    " < {snakemake.input.call} > {snakemake.output} {log}"
)
