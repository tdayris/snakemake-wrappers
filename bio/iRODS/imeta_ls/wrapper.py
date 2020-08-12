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

if "collection" in snakemake.input.keys():
    with open(snakemake.input["collection"], "r") as collection_list:
        for line in collection_list:
            collection = line[:-1]
            shell(
                "imeta ls "  # iRODS command
                " {extra} "  # Extra parameters
                " {collection} "  # Name of the collection
                " {attribute} "  # Name of the attribute to search for
                " >> {snakemake.output[0]} "  # Path to output file
                " {log} "  # Logging behavior
            )

if "name" in snakemake.params.keys():
    for collection_name in snakemake.params["name"]:
        shell(
            "imeta ls "  # iRODS command
            " {extra} "  # Extra parameters
            " {collection_name} "  # Name of the collection
            " {attribute} "  # Name of the attribute to search for
            " >> {snakemake.output[0]} "  # Path to output file
            " {log} "  # Logging behavior
        )
