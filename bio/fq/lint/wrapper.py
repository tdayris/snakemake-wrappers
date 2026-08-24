# coding: utf-8

"""snakemake wrapper for fq lint"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

fq_nb = len(snakemake.input)
if not (1 <= fq_nb <= 2):
    raise ValueError(f"This wrapper expected 1 or 2 fastq files, got {fq_nb}.")

shell("fq lint {extra} {snakemake.input} > {snakemake.output:q} {log}")
