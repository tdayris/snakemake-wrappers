# coding: utf-8

"""Snakemake wrappers for sha256sum"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from shlex import quote

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("sha256sum {snakemake.input} > {snakemake.output:q} {log}")
