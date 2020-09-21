__author__ = "Fabian Kilpert"
__copyright__ = "Copyright 2020, Fabian Kilpert"
__email__ = "fkilpert@gmail.com"
__license__ = "MIT"


from snakemake.shell import shell


log = snakemake.log_fmt_shell(stdout=True, stderr=True)


memory = ""
if "mem_mb" in snakemake.resources.keys():
    memory = "-Xmx{}M".format(snakemake.resources["mem_mb"])


shell(
    "picard BedToIntervalList "
    " {memory} "
    " {snakemake.params} "
    " INPUT={snakemake.input.bed} "
    " SEQUENCE_DICTIONARY={snakemake.input.dict} "
    " OUTPUT={snakemake.output} "
    " {log} "
)
