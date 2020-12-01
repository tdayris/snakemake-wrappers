#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for rust-bio-tools sequence-stats"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

min_threads = 1
extra = snakemake.params.get("extra", "")
sequence = snakemake.input.get("seq", "")

# Adding the -q flag in case user provided a Fastq file
if sequence.endswith(("fq", "fastq", "fq.gz", "fastq.gz")):
    extra += " --fastq "

# Dynamically gunzipping sequence file if needed
if sequence.endswith(".gz"):
    sequence = "<(gunzip -c {})".format(sequence)
    min_threads += 1  # 1 rbt + 1 gunzip = 2 threads


if snakemake.threads < min_threads:
    raise ValueError(
        "Expected at least {} threads, got {}".format(
            min_threads, snakemake.threads
        )
    )


shell("rbt sequence-stats {extra} < {sequence} > {snakemake.output} {log}")
