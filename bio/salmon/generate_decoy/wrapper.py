#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for bcl2fastq"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)
min_threads = 2

genome_input = snakemake.input["genome"]
if genome_input.endswith(".gz"):
    genome_input = "<(gunzip -c {})".format(genome_input)
    min_threads += 1

if snakemake.threads != min_threads:
    raise ValueError(
        "This wrapper requires exactly {} threads".format(min_threads)
    )

shell(
    "grep '^>' {genome_input} | cut -d ' ' -f 1 "
    "> {snakemake.output.decoys} {log}"
)

shell(
    "sed -i.bak -e 's/>//g' {snakemake.output.decoys} {log} & "
    "cat {snakemake.input.transcriptome} {snakemake.input.genome} "
    "> {snakemake.output.gentrome} {log} & wait"
)
