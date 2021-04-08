#!/bin/bash -e

# Authors: Thibault Dayris
# Copyright: Copyright 2021, Thibault Dayris
# Email: thibault.dayris@gustaveroussy.fr
# License: MIT

# This script is used to export variables to a running environment

# Add shortcut to conda environment
declare -x CONDA_ENV_PATH="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/env/"

# Declare snakemake cache directory
declare -x SNAKEMAKE_OUTPUT_CACHE="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/cache"

# Export previously defined variables to current environment
export SNAKEMAKE_OUTPUT_CACHE CONDA_ENV_PATH
