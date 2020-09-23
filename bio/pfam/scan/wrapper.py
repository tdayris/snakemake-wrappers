#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for Pfam scan"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
from os.path import dirname
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

pfam_dir = snakemake.input.pfam_dir

shell(
    " pfam_scan.pl "
    " -fasta {snakemake.input.fasta} "
    " -dir {pfam_dir} "
    " --outfile {snakemake.output[0]} "
    " -cpu {snakemake.threads} "
    " {extra} "
    " {log} "
)
