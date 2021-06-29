"""
Snakemake wrapper for GATK CreateSomaticPanelOfNormals
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import os

from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts


extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

gvcfs = list(map("--variant {}".format, snakemake.input.gvcfs))

bam = ""
if "bam" in snakemake.input.keys():
    bam = list(map("--intervals {}".format, snakemake.input.bam))

intervals = ""
if "intervals" in snakemake.input.keys():
    intervals = list(map("--input {}".format, snakemake.input.intervals))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "gatk --java-options '{java_opts}' CreateSomaticPanelOfNormals {extra} "
    "{gvcfs} {bam} {intervals} "
    "--reference {snakemake.input.ref} "
    "--output {snakemake.output.gvcf} {log}"
)
