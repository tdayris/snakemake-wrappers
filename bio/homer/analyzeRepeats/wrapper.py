#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Snakemake wrapper for analyzeRepeats.pl"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2022, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Defining input type/path
tag_directories = " ".join(snakemake.input["tag"])

data_input = snakemake.params.get("gene_definition", None)
if "gtf" in snakemake.input.keys():
    if data_input is not None:
        raise ValueError(
            "Input should be 'rna', 'repeats', "
            "a path to a GTF file **OR** a peaks file."
        )
    else:
        data_input = snakemake.input["gtf"]

if "peaks" in snakemake.input.keys():
    if data_input is not None:
        raise ValueError(
            "Input should be 'rna', 'repeats', "
            "a path to a GTF file **OR** a peaks file."
        )
    else:
        data_input = snakemake.input["peaks"]


# Defining genome version and other optional parameters
genome = snakemake.params.get("genome", None)
if "fasta" in snakemake.input.keys():
    if genome is not None:
        raise ValueError(
            "Genome should be a string **OR** a path to a "
            "fasta file. Not both."
        )
    else:
        genome = snakemake.input["fasta"]

extra = snakemake.params.get("extra", "")


shell(
    "analyzeRepeats.pl "
    "{data_input} "
    "{genome} "
    "{extra} "
    " -d {tag_directories} "
    " > {snakemake.output[0]} "
    "{log}"
)
