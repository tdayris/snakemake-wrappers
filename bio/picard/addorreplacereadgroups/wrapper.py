__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

memory = ""
if "mem_mb" in snakemake.resources.keys():
    memory = "-Xmx{}M".format(snakemake.resources["mem_mb"])


shell(
    " picard AddOrReplaceReadGroups "
    " {memory} "
    " {snakemake.params} "
    " I={snakemake.input} "
    " O={snakemake.output} "
    " {log}"
)
