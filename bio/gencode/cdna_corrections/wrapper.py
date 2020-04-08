#!/bin/python3.8
# conding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

"""
Corrects transcripts names from the CDNA fasta file from gencode
"""

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)
input_stream =  snakemake.input["fasta"]
output_stream = f"> {snakemake.output['fasta']}"
min_threads = 1

if snakemake.input["fasta"].endswith("gz"):
    input_stream = f"<(gunzip -c {snakemake.input['fasta']})"
    min_threads += 1

if snakemake.output["fasta"].endswith("gz"):
    output_stream = f"| gzip -c > {snakemake.output['fasta']})"
    min_threads += 1

if snakemake.threads < min_threads:
    raise ValueError(
        f"At least {min_threads} are required for gencode corretions"
    )


shell("sed 's/|.*//g' {input_stream} {output_stream} {log}")
