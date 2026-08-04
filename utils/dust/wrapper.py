# coding: utf-8

"""Snakemake wrapper used for dust"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import get_format

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

if get_format(snakemake.output[0]) == "json":
    extra += " --output-json"

shell(
    "dust {extra} --threads {snakemake.threads} "
    "{snakemake.input.path} > {snakemake.output} {log}"
)
