# coding: utf-8

"""Snakemake wrapper for fclones group command"""

from os.path import isfile, islink
from snakemake.shell import shell
from snakemake_wrapper_utils.snakemake import get_format

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra = snakemake.params.get("extra", "")

out_format = get_format(str(snakemake.output[0]))
if out_format in ("csv", "json"):
    extra += f" --format '{out_format}'"

if all(isfile(str(f)) for f in snakemake.input.paths):
    extra += " --depth 0"

if any(str(f).startswith(".") for f in snakemake.input.paths):
    extra += " --hidden"

if any(islink(str(f)) for f in snakemake.input.paths):
    extra += " --symbolic-links"

shell("fclones group {extra} {snakemake.input.paths} --output {snakemake.output[0]:q} {log}")
