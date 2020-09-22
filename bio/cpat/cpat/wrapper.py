#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for CPAT"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# A reference sequence in FASTA format is needed only if
# transcript file was in BED format.
reference = ""
if "reference" in snakemake.input.keys():
    reference = "--ref {}".format(snakemake.input["reference"])

extra = snakemake.params.get("extra", "")

shell(
    " cpat.py "
    " --gene {snakemake.input.transcripts} "
    " {reference} "
    " --hex {snakemake.input.hexamer_table} "
    " --logitModel {snakemake.input.logit_model}"
    " --outfile {snakemake.output[0]}"
    " {extra} "
    " {log} "
)
