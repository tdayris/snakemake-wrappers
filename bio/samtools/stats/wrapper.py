"""Snakemake wrapper for trimming paired-end reads using cutadapt."""

__author__ = "Julian de Ruiter"
__copyright__ = "Copyright 2017, Julian de Ruiter"
__email__ = "julianderuiter@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell


extra = snakemake.params.get("extra", "")
region = snakemake.params.get("region", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

reference = snakemake.input.get("ref", "")
if reference:
    reference = f"-r {reference}"


regions = snakemake.input.get("bed", "")
if regions:
    regions = f"-t {regions}"


shell("samtools stats {extra} {reference} {regions} {snakemake.input.aln} {region} > {snakemake.output} {log}")
