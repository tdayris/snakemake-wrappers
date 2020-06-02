"""Snakemake wrapper for SnpSift annotate"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from os.path import dirname
from snakemake.shell import shell
from snakemake.utils import makedirs

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

extra = snakemake.params.get("extra", "")

incall = snakemake.input.incall
if incall.endswith(".bcf"):
    incall = "< <(bcftools view {})".format(incall)

shell(
    "SnpSift annotate"  # Tool and its subcommand
    " {extra}"  # Extra parameters
    " {snakemake.input.database}"  # Path to annotation vcf file
    " {incall}"  # Path to input vcf file
    " > {snakemake.output.vcf}"  # Path to output vcf file
    " {log}"  # Logging behaviour
)
