#!/bin/bash -e

# Authors: Thibault Dayris
# Copyright: Copyright 2021, Thibault Dayris
# Email: thibault.dayris@gustaveroussy.fr
# License: MIT

# This script is used to export variables to a running environment

function conda_activate () {
  # shellcheck source=/mnt/beegfs/userdata/t_dayris/anaconda3/etc/profile.d/conda.sh
  source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate && conda activate "${1}"
}

function profiles () {
  local PROFILE_NAME="${1}"
  message INFO "Looking for profile named: ${PROFILE_NAME}"
  case "${PROFILE_NAME}" in
    demux) echo "${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/demux";;
    slurm) echo "${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/slurm";;
    clinics) echo "${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/clinics";;
    *) echo "Unknown profile name";;
  esac
}

function iodirectories() {
  CWD="${1}"
  DIRNAME="${2}"

  # Default IO directories
  DIRNAME_PATH="$(readlink -e ${CWD})/${DIRNAME}"
  if [ -L "${DIRNAME_PATH}" ]; then
    if [ -e "${DIRNAME_PATH}" ]; then
      # Case DIRNAME is a symlink
      message INFO "${DIRNAME} symlink directory already available: ${DIRNAME_PATH}"
    else
      message ERROR "${DIRNAME} is a broken symlink: ${DIRNAME_PATH}. This will raise error in the future."
    fi
  elif [ -d "${DIRNAME_PATH}" ]; then
    # case DIRNAME is a dir
    message WARNING "${DIRNAME} directory already available: ${DIRNAME_PATH}. It does not seem to be linked to official BiGR data managment"
  elif [ -d $(readlink -e "${CWD}/../${DIRNAME}") ]; then
    # Case DIRNAME is missing but available in parent dir
    ln -sfrv $(readlink -e "${CWD}/../${DIRNAME}") "${DIRNAME_PATH}"
    message INFO "${DIRNAME} directory available from parent dir and linked to: ${DIRNAME_PATH}"
  else
    # Case DIRNAME never found
    message WARNING "Data input was not found. It will be created by snakemake, but not linked to official BiGR data managment"
  fi
}

# Add shortcut to conda environment, the main environment with resources to execute all pipelines
declare -x CONDA_ENV_PATH="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/bigr_snakemake"

# Declare snakemake cache directory. Used to avoid indexation steps and redundant operations
declare -x SNAKEMAKE_OUTPUT_CACHE="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/snakemake_cache"

# Declare conda cache directory. Used to avoid conda reinstallations
declare -x CONDA_CACHE_PATH="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/conda_cache"

# Export previously defined variables to current environment
export SNAKEMAKE_OUTPUT_CACHE CONDA_ENV_PATH CONDA_CACHE_PATH

mkdir --parents --verbose tmp/shadow

message INFO "sourced: environment.sh"
