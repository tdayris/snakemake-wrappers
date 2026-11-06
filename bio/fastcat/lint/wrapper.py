# coding: utf-8

"""Snakemake wrapper for running quality control with fastcat lint and run sdust"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("fastcat lint {extra} {snakemake.input[0]:q} > {snakemake.output[0]:q} {log}")
