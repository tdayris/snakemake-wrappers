__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts


samtools_opts = get_samtools_opts(snakemake)
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


tmp_dir = snakemake.params.get("tmp_dir", "")
if tmp_dir:
    prefix = os.path.join(tmp_dir, os.path.basename(out_name))
else:
    prefix = out_name

if not os.path.exists(tmp_dir):
    os.makedirs(name=tmp_dir, exist_ok=True)

# Samtools takes additional threads through its option -@
# One thread for samtools
# Other threads are *additional* threads passed to the argument -@
threads = "" if snakemake.threads <= 1 else " -@ {} ".format(snakemake.threads - 1)

# Memory per thread
extra = snakemake.params[0]
if "mem_mb" in snakemake.resources.keys() and "-m" not in extra:
    extra += " -m {}M".format(int(snakemake.resources["mem_mb"] / snakemake.threads))

shell(
    "samtools sort {snakemake.params.extra} {threads} -o {snakemake.output[0]} "
    "-T {prefix} {snakemake.input[0]} {log}"
)
