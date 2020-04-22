"""Snakemake wrapper for SnpSift dbNSFP"""

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

db = snakemake.input.get("dbNSFP", "")
if db != "":
    db = "-db {}".format(db)

shell(
    "SnpSift dbnsfp"  # Tool and its subcommand
    " {extra}"  # Extra parameters
    " {db}"  # Path to annotation vcf file
    " -v {snakemake.input.vcf}"  # Path to input vcf file
    " > {snakemake.output.vcf}"  # Path to output vcf file
    " {log}"  # Logging behaviour
)
