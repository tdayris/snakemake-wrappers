#!/bin/python3.8

# This wrapper downloads and install optional EaCoN packages

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.utils import makedirs
from snakemake.shell import shell
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

for package in snakemake.input.get("r_packages", []):
     shell(
        'R  --vanilla CMD INSTALL --clean --data-compress=gzip {package} {log}'
    )

for package in snakemake.input.get("git_packages", []):
    shell(
        'R --vanilla -e \'options(unzip = "internal"); '
        'Sys.setenv(TAR = \"$(which tar)\"); '
        'devtools::install("{package}", force=TRUE);\' {log}'
    )
