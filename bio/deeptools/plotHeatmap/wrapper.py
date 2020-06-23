#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""This is the Snakemake Wrapper for deeptools plotHeatmap"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

extra = snakemake.params.get("extra", "")

if "matrix" in snakemake.output.keys():
    extra += " --outFileNameMatrix {} ".format(snakemake.output.matrix)
if "regions" in snakemake.output.keys():
    extra += " --outFileSortedRegions {} ".format(snakemake.output.regions)

shell(
    "plotHeatmap "
    " --matrixFile {snakemake.input.matrix} "
    " --outFileName {snakemake.output.plot} "
    " {extra} "
    " {log} "
)
