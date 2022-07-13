"""Snakemake wrapper for Varscan2 mpileup2snp"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Dayris Thibault"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os.path as op
from snakemake.shell import shell
from snakemake.utils import makedirs
from snakemake_wrapper_utils.java import get_java_opts

# Gathering extra parameters and logging behaviour
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

# In case input files are gzipped mpileup files,
# they are being unzipped and piped
# In that case, it is recommended to use at least 2 threads:
# - One for unzipping with zcat
# - One for running varscan
pileup = (
    " cat {} ".format(snakemake.input["pileup"])
    if not snakemake.input["pileup"].endswith("gz")
    else " zcat {} ".format(snakemake.input["pileup"])
)

if str(snakemake.output[0]).endswith(("vcf", "vcf.gz")):
    extra += " --output-vcf 1 "

if "sample_list" in snakemake.input.keys():
    extra += " --vcf-sample-list {} ".format(snakemake.input['sample_list'])

# Building output directories
makedirs(op.dirname(snakemake.output[0]))

shell(
    "varscan mpileup2cns "  # Tool and its subprocess
    "<( {pileup} ) "
    "{java_opts} {extra} "  # Extra parameters
    "> {snakemake.output[0]} "  # Path to vcf file
    "{log}"  # Logging behaviour
)
