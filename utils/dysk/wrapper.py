# coding: utf-8

"""Snakemake wrapper for dysk"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import get_format

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

outfile_ext = get_format(snakemake.output[0])
if outfile_ext == "json":
    extra += f" --json"
elif outfile_ext == "csv":
    extra += f" --csv"
elif outfile_ext == "tsv":
    extra += f" --csv --csv-separator $'\t'"

shell("dysk {extra} {snakemake.input} > {snakemake.output} {log}")
