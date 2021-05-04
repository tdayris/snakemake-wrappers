__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params
java_opts = get_java_opts(snakemake)

shell(
    "picard CollectAlignmentSummaryMetrics INPUT={snakemake.input.bam} OUTPUT={snakemake.output[0]} "
    "REFERENCE_SEQUENCE={snakemake.input.ref} {log}"
)
