#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper for bash copy"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os.path as op
import yaml

from snakemake.shell import shell

# No node can access a cold storage
# these files must be copied. However,
# any file else where should be symlinked!
cold_storage = None
if isinstance(snakemake.params.cold_storage, list):
    cold_storage = snakemake.params.cold_storage
elif isinstance(snakemake.params.cold_storage, str):
    if op.exists(snakemake.params.cold_storage):
        with open(snakemake.params.cold_storage, "r") as cold_list:
            cold_storage = yaml.load(cold_list)
    else:
        cold_storage = [snakemake.params.cold_storage]
else:
    cold_storage = [str(snakemake.params.cold_storage)]

# Prepare logging
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
extra_cp = snakemake.params.get("extra", "-v")
extra_iget = snakemake.params.get("extra_irods", "-vK")

output_directory = op.realpath(op.dirname(snakemake.output[0]))

shell_cmd = f"mkdir --parents {output_directory} && "
for f in snakemake.input:
    if f.startswith("/odin/kdi/dataset"):
        shell_cmd += (
            "iget "                                      # iRODS copy
            " {extra_iget} "                             # Extra parameters
            f"{op.realpath(f)} "                         # Path to input file
            f"{op.realpath(snakemake.output[0])} "       # Path to output file
            "{log}; "                                    # Logging
        )
    else:
        symlink = any(f.startswith(cs) for cs in cold_storage)
        shell_cmd += (
            "cp "                                        # Tool
            "{extra_cp} "                                # Optionnal parameters
            f"{'' if symlink else '--symbolic-link'} "   # Symlink if needed
            f"{op.realpath(f)} "                         # Path to input file
            f"{op.realpath(snakemake.output[0])} "       # Path to output file
            "{log}; "                                    # Logging
        )

shell(shell_cmd)
