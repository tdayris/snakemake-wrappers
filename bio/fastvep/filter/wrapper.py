# coding: utf-8

"""Snakemake wrapper for fastvep external annotation source build"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from shlex import quote

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "fastvep filter {extra} --input {snakemake.input.vcf:q} "
    "--output {snakemake.output.vcf:q} {log}"
)
