# coding: utf-8

"""Snakemake wrapper for fastvep external annotation source build"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2026, Thibault Dayris"
__license__ = "MIT"

from os.path import commonprefix
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

output_extensions = ("osa", "osa2", "osi", "oga")
outfile_prefix = str(snakemake.output[0])
for ext in output_extensions:
    if outfile_prefix.endswith(ext):
        outfile_prefix = outfile_prefix[: -len(ext) + 1]
        if ext in outfile_prefix[:2]:
            extra += f" --format '{ext}'"
        break

shell(
    "fastvep sa-build {extra} --input {snakemake.input[0]:q} "
    "--output {outfile_prefix:q} {log}"
)
