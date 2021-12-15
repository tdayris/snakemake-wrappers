__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
import tempfile

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

mapin = snakemake.input["mapping"]
mapout = snakemake.output["mapping"]
tempdir = tempfile.mkdtemp()
if "tmpdir" in snakemake.resources.keys():
    tempdir = snakemake.resources["tmpdir"]

memory = "--memory-limit "
if "mem_mb" in snakemake.resources.keys():
    memory += str(int(snakemake.resources["mem_mb"] / snakemake.threads)) + "MB"
else:
    memory += "2GB"


shell(
    "sambamba sort "
    "{snakemake.params} "
    "{memory} "
    "--tmpdir {tempdir} "
    "--nthreads {snakemake.threads} "
    "--out {mapout} "
    "{mapin} "
    "{log}"
)
