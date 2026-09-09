# coding: utf-8

"""Snakemake wrappers used to run bam_annotate"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from shlex import quote

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

genes = snakemake.input.get("genes")
if genes:
    extra += f" -g {quote(genes)}"

transcripts = snakemake.input.get("transcripts")
if transcripts:
    extra += f" -t {quote(transcripts)}"

introns = snakemake.input.get("introns")
if introns:
    extra += f" -i {quote(introns)}"

exons = snakemake.input.get("exons")
if exons:
    extra += f" -e {quote(exons)}"

shell(
    "bam_annotate.sh -b {snakemake.input.aln:q} {extra} "
    "> {snakemake.output[0]:q} {log}"
)
