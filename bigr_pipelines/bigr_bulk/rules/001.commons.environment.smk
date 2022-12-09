import os

from snakemake.utils import min_version

min_version("7.18")

envvars:
    "SNAKEMAKE_OUTPUT_CACHE"