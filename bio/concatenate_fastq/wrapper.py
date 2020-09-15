#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

"""Snakemake wrapper fastq concatenation"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2019, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

if snakemake.threads != 2:
    raise ValueError(
        "Exactly two threads are required for this wrapper. "
        "Only {} were given.".format(str(snakemake.threads))
    )

if len(snakemake.input) > 1:
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)
    shell(
        " zcat "  # Concatenate gzipped files
        " {snakemake.input} "  # List of input files
        " | "
        " gzip -c "  # Gzip output file
        " > {snakemake.output} "  # Path to output file
        " {log} "  # Logging
    )
else:
    log = snakemake.log_fmt_shell(stdout=True, stderr=True)
    shell(
        " cp -v "  # Copy if not merge needed
        " {input} "  # Path to input file
        " {output} "  # Path to output file
        " {log} "  # Logging
    )
