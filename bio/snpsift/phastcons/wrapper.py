"""Snakemake wrapper for SnpSift dbNSFP"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

min_threads = 1

# Uncompression shall be done on user request
incall = snakemake.input["call"]
if incall.endswith("bcf"):
    min_threads += 1
    incall = "< <(bcftools view {})".format(incall)
elif incall.endswith("gz"):
    min_threads += 1
    incall = "< <(gunzip -c {})".format(incall)

# Compression shall be done according to user-defined output
outcall = snakemake.output["call"]
if outcall.endswith("gz"):
    min_threads += 1
    outcall = "| bcftools view --output-type z > {}".format(outcall)
elif outcall.endswith("bcf"):
    min_threads += 1
    outcall = "| bcftools view --output-type b > {}".format(outcall)
else:
    outcall = "> {}".format(outcall)

# Each (un)compression raises the thread number
if snakemake.threads < min_threads:
    raise ValueError(
        "At least {} threads required, {} provided".format(
            min_threads, snakemake.threads
        )
    )


shell(
    "SnpSift phastCons"  # Tool and its subcommand
    " {java_opts} {extra}"  # Extra parameters
    " {snakemake.input.phastcons}"  # Path to annotation directory
    " {incall}"  # Path to input vcf file
    " {outcall}"  # Path to output vcf file
    " {log}"  # Logging behaviour
)
