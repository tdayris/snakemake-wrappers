__author__ = "Jan Forster"
__copyright__ = "Copyright 2021, Jan Forster"
__email__ = "j.forster@dkfz.de"
__license__ = "MIT"


import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

mapping = snakemake.input["mapping"]
if mapping.endswith(("fa", "fasta")) and "fasta" not in extra:
        extra += " --fasta-input"

regions = snakemake.input.get("regions", "")
if regions != "" and "regions" not in extra:
    extra += " --regions {}".format(regions)

shell("sambamba slice {mapping} {extra} > {snakemake.output[0]} {log}")
