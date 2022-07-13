"""Snakemake wrapper for Salmon Index."""

__author__ = "Tessa Pierce"
__copyright__ = "Copyright 2018, Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

if "decoys" in snakemake.input.keys():
    extra += " --decoys {} ".format(snakemake.input["decoys"])

shell(
    "salmon index --transcripts {snakemake.input.sequences} "
    "--index {snakemake.output.index} "
    "--threads {snakemake.threads} {extra} {log}"
)
