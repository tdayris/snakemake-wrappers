# coding: utf-8

"""Snakemake wrapper for fclones group command"""

from os.path import isfile, islink
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

shell("grex {extra} --file {snakemake.input.test_cases:q} > {snakemake.output:q} {log}")
