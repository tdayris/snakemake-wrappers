#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for CPAT make_legit_model.py"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")
prefix = snakemake.output["rdata"][:-len(".logit.RData")]

reference = ""
if "reference" in snakemake.input.keys():
    reference = "--ref {snakemake.input.reference}"

shell(
    " make_logitModel.py "
    " --cgene {snakemake.input.coding} "
    " --ngene {snakemake.input.noncoding} "
    " --hex {snakemake.input.hexamer_table} "
    " --outfile {prefix} "
    " {reference} "
    " {extra} "
    " {log} "
)
