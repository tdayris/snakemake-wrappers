#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper for iRODS imeta qu"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell

# Prepare logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra = snakemake.params.get("extra", "")

kwords_list = [
    f'"{k}" {"like" if "%" in v else "="} "{v}"'
    for k, v in snakemake.params["kwords"].items()
]
kwords_command = " and ".join(kwords_list)

shell(
    "imeta qu"  # iRODS command
    " {extra} "  # Extra parameters
    " {kwords_command} "  # Kwords list
    " > {snakemake.output[0]} "  # Path to output file
    " {log} "  # Logging behavior
)
