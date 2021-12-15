#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Sansa wrapper for Snakemake"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

gtf = ""
if "gtf" in snakemake.input.keys():
    gtf = "--gtf {}".format(snakemake.input["gtf"])

db = ""
if "db" in snakemake.input.keys():
    db = "--db {}".format(snakemake.input["db"])

contained = ""
if "contained" in snakemake.output.keys():
    contained = "--contained {}".format(snakemake.output["contained"])

#anno = "--anno {}".format(snakemake.output["vcf"])
extra = snakemake.params.get("extra", "")

shell(
    "sansa annotate "
    "{gtf} "
    "{db} "
    #"{anno} "
    "{contained} "
    "{extra} "
    "{snakemake.input.vcf} "
    " > {snakemake.output.vcf} "
    "{log}"
)
