__author__ = "Jan Forster"
__copyright__ = "Copyright 2021, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"


import os

from snakemake.shell import shell

in_file = snakemake.input["mapping"]

regions = snakemake.input.get("regions", "")
extra = snakemake.params.get("extra", "")

if regions != "":
    extra = "--regions {}".format(regions)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if in_file.endswith(".sam") and ("-S" not in extra or "--sam-input" not in extra):
    extra += " --sam-input"

shell(
    "sambamba view {extra} --nthreads {snakemake.threads} "
    "{in_file} --output-filename {snakemake.output[0]} {log}"
)
