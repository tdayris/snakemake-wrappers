#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for mixcr assemble"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")


if "report" in snakemake.output.keys():
    extra += " --report {} ".format(snakemake.output["report"])


shell(
    "mixcr assemble {extra} --threads {snakemake.threads} "
    "{snakemake.input} {snakemake.output.clones} {log}"
)
