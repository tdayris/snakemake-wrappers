#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper for iRODS imeta ls"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell

# Prepare logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)
extra = snakemake.params.get("extra", "")
attribute = snakemake.params.get("attribute", "")

if isinstance(snakemake.params["name"], str):
    shell(
        "imeta ls "  # iRODS command
        " {extra} "  # Extra parameters
        " {snakemake.params['name']} "  # Name of the collection
        " {attribute} "  # Name of the attribute to search for
        " > {snakemake.output[0]} "  # Path to output file
        " {log} "  # Logging behavior
    )
else:
    for idx, name in enumerate(snakemake.params['name']):
        overwrite = ">" if idx == 0 else ">>"
        shell(
            "imeta ls "  # iRODS command
            " {extra} "  # Extra parameters
            " {name} "  # Name of the collection
            " {attribute} "  # Name of the attribute to search for
            " {overwrite} {snakemake.output[0]} "  # Path to output file
            " {log} "  # Logging behavior
        )
