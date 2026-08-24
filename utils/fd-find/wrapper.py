# coding: utf-8

"""Snakemake wrapper for fd-find"""

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("fd {extra} --threads {snakemake.threads} {snakemake.input:q} > {snakemake.output:q} {log}")
