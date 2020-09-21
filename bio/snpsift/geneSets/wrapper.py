"""Snakemake wrapper for SnpSift geneSets"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os.path import dirname
from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")
makedirs(dirname(snakemake.output.vcf))


if "mem_mb" in snakemake.resources.keys():
    extra += "-Xmx{}M".format(snakemake.resources["mem_mb"])

shell(
    "SnpSift geneSets"  # Tool and its subcommand
    " {extra}"  # Extra parameters
    " {snakemake.input.gmt}"  # Path to annotation vcf file
    " -v {snakemake.input.vcf}"  # Path to input vcf file
    " > {snakemake.output.vcf}"  # Path to output vcf file
    " {log}"  # Logging behaviour
)
