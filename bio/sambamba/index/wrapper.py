__author__ = "Jan Forster"
__copyright__ = "Copyright 2021, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"


import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if snakemake.input[0].endswith(("fa", "fasta")):
    extra += "--fasta-input"

shell(
    "sambamba index {snakemake.params.extra} --nthreads {snakemake.threads} "
    "{snakemake.input[0]} {snakemake.output[0]} "
    "{log}"
)
