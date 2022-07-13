__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
import tempfile

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "sambamba sort {snakemake.params} --nthreads {snakemake.threads} "
        "--out {snakemake.output[0]} {snakemake.input[0]} "
        "--tmpdir {tmpdir}"
        "{log}"
    )
