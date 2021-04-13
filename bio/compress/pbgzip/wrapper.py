__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2021, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell(
    "pbgzip -c {snakemake.params} "
    "{snakemake.input[0]} "
    "> {snakemake.output[0]} {log}"
)
