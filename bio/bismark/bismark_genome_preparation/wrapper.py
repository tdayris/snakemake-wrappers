"""Snakemake wrapper for Bismark indexes preparing using bismark_genome_preparation."""
# https://github.com/FelixKrueger/Bismark/blob/master/bismark_genome_preparation

__author__ = "Roman Chernyatchik"
__copyright__ = "Copyright (c) 2019 JetBrains"
__email__ = "roman.chernyatchik@jetbrains.com"
__license__ = "MIT"


from os import path
from snakemake.shell import shell

input_dir = path.dirname(snakemake.input[0])

params_extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# From documentation:
# Remember that the indexing is run twice in parallel already
# so e.g. '--parallel 4' will use 8 threads in total.
parallel = f"--parallel {int(snakemake.threads // 2)}"

shell(
    "bismark_genome_preparation "
    "{parallel} "
    "--verbose "
    "--bowtie2 "
    "{params_extra} "
    "{input_dir} "
    "{log}")
