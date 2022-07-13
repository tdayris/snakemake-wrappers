"""Snakemake wrapper for subsampling reads from paired FASTQ files using seqtk."""

__author__ = "Fabian Kilpert"
__copyright__ = "Copyright 2020, Fabian Kilpert"
__email__ = "fkilpert@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell


log = snakemake.log_fmt_shell()
pigz_threads = snakemake.threads - 1
if snakemake.threads <= 1:
    raise ValueError("This wrapper requires at least two threads")

shell(
    "( "
    "seqtk sample "
    "-s {snakemake.params.seed} "
    "{snakemake.input.f1} "
    "{snakemake.params.n} "
    "| pigz -p {pigz_threads} "
    "> {snakemake.output.f1} "
    "&& "
    "seqtk sample "
    "-s {snakemake.params.seed} "
    "{snakemake.input.f2} "
    "{snakemake.params.n} "
    "| pigz -p {pigz_threads} "
    "> {snakemake.output.f2} "
    ") {log} "
)
