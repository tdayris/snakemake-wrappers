# coding: utf-8

"""Snakemake wrapper for rbt vcf-fix-iupac-alleles"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake_wrapper_utils.bcftools import get_bcftools_opts
from snakemake_wrapper_utils.snakemake import get_format
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

write_cmd = f"> {snakemake.output[0]} "
if get_format(snakemake.output[0]) != "bcf":
    bcftools_write_opts = get_bcftools_opts(snakemake, parse_memory=False)
    write_cmd = f"| bcftools view {bcftools_write_opts} "

if snakemake.threads < 2:
    raise ValueError("This wrapper expects at least two threads.")

shell("(rbt vcf-fix-iupac-alleles < {snakemake.input[0]} {write_cmd} ) {log}")
