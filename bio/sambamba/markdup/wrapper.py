__author__ = "Jan Forster"
__copyright__ = "Copyright 2021, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"


import os
import tempfile
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

tmpdir = tempfile.mkdtemp()
if "tmpdir" in snakemake.resources.keys():
    tmpdir = snakemake.resources["tmpdir"]

shell(
    "sambamba markdup {snakemake.params.extra} --nthreads {snakemake.threads} "
    "{snakemake.input[0]} {snakemake.output[0]} --tmpdir {tmpdir}"
    "{log}"
)
