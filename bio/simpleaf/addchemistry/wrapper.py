__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thiault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from os.path import basename

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

shell(
    "simpleaf add-chemistry "
    "--name {snakemake.params.name} "
    "--geometry {snakemake.params.geometry} "
    "{log} "
)

if basename(snakemake.output[0]) != "custom_chemistries.json":
    shell("mv --verbose custom_chemistries.json {snakemake.output[0]} {log}")