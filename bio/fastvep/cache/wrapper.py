# coding: utf-8

"""Snakemake wrapper for fastvep external annotation source build"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from shlex import quote

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

synonyms = snakemake.input.get("synonyms")
if synonyms:
    extra += f" --synonyms {quote(synonyms)}"

shell(
    "fastvep cache {extra} "
    "--gff3 {snakemake.input.gff:q} "
    "--fasta {snakemake.input.fasta:q} "
    "--output {snakemake.output[0]:q} "
    "{log}"
)
