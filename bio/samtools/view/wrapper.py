__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os

from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts

samtools_opts = get_samtools_opts(snakemake)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

if "aln" in snakemake.input.keys():
    input_file = snakemake.input["aln"]
else:
    input_file = snakemake.input[0]

shell("samtools view {samtools_opts} {extra} {input_file} {log}")
