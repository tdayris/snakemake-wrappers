#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper designed to extract """

__author__ = "Thibault Dayris"
__copyright__ = "Fish_n_CHIP 2019"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "GPLV3"

from snakemake.shell import shell

# Prepare logging
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

if snakemake.threads != 3:
    raise ValueError(
        "The iRODS' extract-collections wrapper requires exactly 3 threads."
    )

shell(
    'grep -P "^(collection|dataObj):"'  # Grep line with interesting content
    ' {snakemake.input[0]} | '  # Path to input file
    ' cut -f2 -d " " | '  # Remove line identifier
    ' paste -d "/" - - '  # Paste collection and dataset names side by side
    ' > {snakemake.output[0]} '  # Path to output file
    ' {log} '  # Logging
)
