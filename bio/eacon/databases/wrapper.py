#!/bin/python3.8

# This wrapper downloads and install optional EaCoN packages

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

if "databases" in snakemake.output.keys():
    shell("wget -O {snakemake.output.databases}.tar.bz2 https://nextcloud.gustaveroussy.fr/s/jZkw6cCsAKMGXsz/download {log}")
    shell("tar -xvjf {snakemake.output.databases}.tar.bz2 {log}")

if "grd" in snakemake.output.keys():
    shell("wget -O {snakemake.output.grd} https://raw.githubusercontent.com/tdayris/cel-cnv-eacon/master/scripts/grd {log}")
    shell("chmod u+x {snakemake.output.grd} {log}")
