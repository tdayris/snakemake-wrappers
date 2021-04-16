#!/bin/bash -e

# Authors: Thibault Dayris
# Copyright: Copyright 2021, Thibault Dayris
# Email: thibault.dayris@gustaveroussy.fr
# License: MIT

# This script is used to export variables to a running environment

function conda_activate () {
  source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate && conda activate "${1}"
}

function profiles () {
  local PROFILE_NAME="${1}"
  message INFO "Looking for profile named: ${PROFILE_NAME}"
  case PROFILE_NAME in
    demux) echo "${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/demux";;
    slurm) echo "${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/slurm";;
    clinics) echo "${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/clinics";;
    *) echo "Unknown profile name";;
  esac
}

# Add shortcut to conda environment, the main environment with resources to execute all pipelines
declare -x CONDA_ENV_PATH="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/env/"

# Declare snakemake cache directory. Used to avoid indexation steps and redundant operations
declare -x SNAKEMAKE_OUTPUT_CACHE="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/cache"

# Export previously defined variables to current environment
export SNAKEMAKE_OUTPUT_CACHE CONDA_ENV_PATH

message INFO "sourced: environment.sh"
