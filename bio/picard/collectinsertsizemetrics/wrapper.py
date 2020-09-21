__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "johannes.koester@protonmail.com"
__license__ = "MIT"


from snakemake.shell import shell


log = snakemake.log_fmt_shell()


memory = ""
if "mem_mb" is snakemake.resources.keys():
    memory = "-Xmx{}M".format(snakemake.resources["mem_mb"])


shell(
    "picard CollectInsertSizeMetrics {memory} {snakemake.params} "
    "INPUT={snakemake.input} OUTPUT={snakemake.output.txt} "
    "HISTOGRAM_FILE={snakemake.output.pdf} {log}"
)
