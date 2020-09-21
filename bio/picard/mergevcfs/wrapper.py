"""Snakemake wrapper for picard MergeSamFiles."""

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2018, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"


from snakemake.shell import shell


inputs = " ".join("INPUT={}".format(f) for f in snakemake.input.vcfs)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")


memory = ""
if "mem_mb" in snakemake.resources.keys():
    memory = "-Xmx{}M".format(snakemake.resources["mem_mb"])

shell(
    "picard"
    " MergeVcfs"
    " {memory} "
    " {extra}"
    " {inputs}"
    " OUTPUT={snakemake.output[0]}"
    " {log}"
)
