# coding: utf-8

"""Snakemake wrappers used to run fastq_info"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2025, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell
from shlex import quote

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

input_fastq = [quote(f) for f in snakemake.input]
if len(input_fastq) == 2:
    extra += " -s"
elif 0 <= len(input_fastq) > 2:
    raise ValueError("This wrapper expected either 1 or 2 fastq files.")

shell("fastq_info {extra} {input_fastq} > {snakemake.output[0]:q} {log}")
