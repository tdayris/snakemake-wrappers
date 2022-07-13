"""Snakemake wrapper for picard CleanSam"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2016, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

shell(
    "picard CleanSam "  # Tool and its subcommand
    "{java_opts} "  # Automatic java option
    "{extra} "  # User defined parmeters
    "INPUT={snakemake.input} "  # Input file
    "OUTPUT={snakemake.output.bam} "  # Output bam
    "{log}"  # Logging
)
