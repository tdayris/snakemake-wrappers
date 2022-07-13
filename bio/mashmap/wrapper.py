#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for MashMap"""

from snakemake.shell import shell

## Extract arguments
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "mashmap --ref {snakemake.input.fasta} "
    "--query {snakemake.input.query} "
    "--threads {snakemake.threads} "
    "--output {snakemake.output} "
    "{extra} "
    "{log}"
)
