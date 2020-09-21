"""Snakemake wrapper for picard RevertSam."""

__author__ = "Patrik Smeds"
__copyright__ = "Copyright 2018, Patrik Smeds"
__email__ = "patrik.smeds@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)


memory = ""
if "mem_mb" is snakemake.resources.keys():
    memory = "-Xmx{}M".format(snakemake.resources["mem_mb"])

shell(
    "picard"
    " RevertSam"
    " {memory} "
    " {extra}"
    " INPUT={snakemake.input[0]}"
    " OUTPUT={snakemake.output[0]}"
    " {log}"
)
