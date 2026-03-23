"""Snakemake wrapper for SnpSift varType"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from snakemake_wrapper_utils.bcftools import get_bcftools_opts

java_opts = get_java_opts(snakemake)
bcftools_opts = get_bcftools_opts(snakemake, parse_ref=False, parse_memory=False)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")
min_threads = 1


# Uncompression shall be done according to user-defined input
incall = snakemake.input[0]
if snakemake.input[0].endswith(("gz", "bcf")):
    min_threads += 1
    incall = "< <(bcftools view {})".format(incall)

# Compression shall be done according to user-defined output
outcall = snakemake.output[0]
if snakemake.output[0].endswith(("gz", "bcf")):
    min_threads += 1
    outcall = "| bcftools view {} > {}".format(bcftools_opts, outcall)
else:
    outcall = "> {}".format(outcall)

# Each (un)compression step raises the threads requirements
if snakemake.threads < min_threads:
    raise ValueError(
        "At least {} threads required, {} provided".format(
            min_threads, snakemake.threads
        )
    )

shell(
    "( SnpSift varType"  # Tool and its subcommand
    " {java_opts} {extra}"  # Extra parameters
    " {incall} {outcall} )"  # Path to input/output vcf files
    " {log}"  # Logging behaviour
)
