"""Snakemake wrapper for SnpSift varType"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

java_opts = get_java_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")


if "mem_mb" in snakemake.resources.keys():
    extra += "-Xmx{}M".format(snakemake.resources["mem_mb"])

shell(
    "SnpSift varType"  # Tool and its subcommand
    " {java_opts} {extra}"  # Extra parameters
    " {snakemake.input.vcf}"  # Path to input vcf file
    " > {snakemake.output.vcf}"  # Path to output vcf file
    " {log}"  # Logging behaviour
)
