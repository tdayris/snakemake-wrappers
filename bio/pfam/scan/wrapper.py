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

pfam_dir = dirname(snakemake.input.pfam_hmm[0])

exts = ["hmm", "hmm.dat", "hmm.h3f", "hmm.h3i", "hmm.h3m", "hmm.h3p"]

shell(
    " pfam_scan.pl "
    " -fasta {snakemake.input.fasta} "
    " -dir {pfam_dir} "
    " --outfile {snakemake.output[0]} "
    " -cpu {snakemake.threads} "
    " {extra} "
    " {log} "
)