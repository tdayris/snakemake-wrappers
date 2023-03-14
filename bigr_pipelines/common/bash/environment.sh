#!/bin/bash -e

# Authors: Thibault Dayris
# Copyright: Copyright 2021, Thibault Dayris
# Email: thibault.dayris@gustaveroussy.fr
# License: MIT

# This script is used to export variables to a running environment

function conda_activate () {
  CMD="source \"/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_conda/anaconda3/etc/profile.d/conda.sh\""
  message CMD "${CMD}"
  eval ${CMD}
  CMD="conda activate"
  message CMD "${CMD}"
  eval ${CMD}
  CMD="conda activate \"${1}\""
  message CMD "${CMD}"
  eval ${CMD}
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
  DIRNAME_PATH="${CWD}/${DIRNAME}"
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
  elif [ -d "${CWD}/../${DIRNAME}" ]; then
    # Case DIRNAME is missing but available in parent dir
    CMD="ln -sfrv \"${CWD}/../${DIRNAME}\" \"${DIRNAME_PATH}\""
    message CMD "${CMD}"
    eval ${CMD}
    message INFO "${DIRNAME} directory available from parent dir and linked to: ${DIRNAME_PATH}"
  else
    # Case DIRNAME never found
    message WARNING "Data input was not found. It will be created, but not linked to official BiGR data managment"
    CMD="mkdir --parents --verbose data_{input,output}"
    message CMD "${CMD}"
    eval ${CMD}
  fi
}

# Add shortcut to conda environment, the main environment with resources to execute all pipelines
declare -x CONDA_ENV_PATH="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/bigr_snakemake"

# Declare snakemake cache directory. Used to avoid indexation steps and redundant operations
declare -x SNAKEMAKE_OUTPUT_CACHE="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/snakemake_cache"

# Declare conda cache directory. Used to avoid conda reinstallations
declare -x CONDA_CACHE_PATH="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/conda_cache"

declare -x SHARED_SINGULARITY_PATH="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/singularity/"
declare -x SHARED_CONDA_INDSTALL="/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_install/"

# Export previously defined variables to current environment
export SNAKEMAKE_OUTPUT_CACHE CONDA_ENV_PATH CONDA_CACHE_PATH

CMD="mkdir --parents --verbose tmp/shadow"
message CMD "${CMD}"
eval ${CMD}

# Overloading TMP directories
# Our main temporary directory.
if [ -z "${BIGR_DEFAULT_TMP:-}" ]; then
  declare -x BIGR_DEFAULT_TMP
  BIGR_DEFAULT_TMP="/mnt/beegfs/userdata/${USER}/tmp"
  export BIGR_DEFAULT_TMP
elif [ "${BIGR_DEFAULT_TMP:-}" == "/tmp" ]; then
  BIGR_DEFAULT_TMP="/mnt/beegfs/userdata/${USER}/tmp"
  export BIGR_DEFAULT_TMP
fi

if [ ! -d "${BIGR_DEFAULT_TMP}" ]; then
  CMD="mkdir --parent --verbose ${BIGR_DEFAULT_TMP}"
  message CMD "${CMD}"
  eval ${CMD}
fi

# Used in many bash / Python scripts
if [ -z "${TMP:-}" ]; then
  message WARNING "TMP environment variable was not set. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "TMP now points to ${BIGR_DEFAULT_TMP}"
  declare -x TMP
  TMP="${BIGR_DEFAULT_TMP}"
  export TMP
elif [ "${TMP:-}" == "/tmp" ]; then
  message WARNING "TMP currently points to '/tmp'. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "This value is now changed to ${BIGR_DEFAULT_TMP}"
  TMP="${BIGR_DEFAULT_TMP}"
  export TMP
else
  message INFO "TMP -> ${TMP}"
fi

# Used in some bash / R / perl / Python scripts
if [ -z "${TEMP:-}" ]; then
  message WARNING "TEMP environment variable was not set. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "TEMP now points to ${BIGR_DEFAULT_TMP}"
  declare -x TEMP
  TEMP="${BIGR_DEFAULT_TMP}"
  export TEMP
elif [ "${TEMP:-}" == "/tmp" ]; then
  message WARNING "TEMP currently points to '/tmp'. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "TEMP now points to ${BIGR_DEFAULT_TMP}"
  TEMP="${BIGR_DEFAULT_TMP}"
  export TEMP
else
  message INFO "TEMP -> ${TEMP}"
fi

# Used in some bash / R / perl / Python scripts
if [ -z "${TMPDIR:-}" ]; then
  message WARNING "TMPDIR environment variable was not set. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "TMPDIR now points to ${BIGR_DEFAULT_TMP}"
  declare -x TMPDIR
  TMPDIR="${BIGR_DEFAULT_TMP}"
  export TMPDIR
elif [ "${TMPDIR:-}" == "/tmp" ]; then
  message WARNING "TMPDIR currently points to '/tmp'. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "TMPDIR now points to ${BIGR_DEFAULT_TMP}"
  TMPDIR="${BIGR_DEFAULT_TMP}"
  export TMPDIR
else
  message INFO "TMPDIR -> ${TMPDIR}"
fi

# Used in some bash / R / perl scripts
if [ -z "${TEMPDIR:-}" ]; then
  message WARNING "TEMPDIR environment variable was not set. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "TEMPDIR now points to ${BIGR_DEFAULT_TMP}"
  declare -x TEMPDIR
  TEMPDIR="${BIGR_DEFAULT_TMP}"
  export TEMPDIR
elif [ "${TEMPDIR:-}" == "/tmp" ]; then
  message WARNING "TEMPDIR currently points to '/tmp'. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "TEMPDIR now points to ${BIGR_DEFAULT_TMP}"
  TEMPDIR="${BIGR_DEFAULT_TMP}"
  export TEMPDIR
else
  message INFO "TEMPDIR -> ${TEMPDIR}"
fi

# Used in nextflow scripts
if [ -z "${NXF_TEMP:-}" ]; then
  message WARNING "NXF_TEMP environment variable was not set. This can lead to NextFlow errors due to lack of space in /tmp"
  message WARNING "NXF_TEMP now points to ${BIGR_DEFAULT_TMP}"
  declare -x NXF_TEMP
  NXF_TEMP="${BIGR_DEFAULT_TMP}"
  export NXF_TEMP
elif [ "${NXF_TEMP:-}" == "/tmp" ]; then
  message WARNING "NXF_TEMP currently points to '/tmp'. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "NXF_TEMP now points to ${BIGR_DEFAULT_TMP}"
  NXF_TEMP="${BIGR_DEFAULT_TMP}"
  export NXF_TEMP
else
  message INFO "NXF_TEMP -> ${NXF_TEMP}"
fi

# Used in nextflow / java scripts
if [ -z "${_JAVA_OPTIONS:-}" ]; then
  message WARNING "_JAVA_OPTIONS environment variable was not set. This can lead to Java errors due to lack of space in /tmp"
  message WARNING "_JAVA_OPTIONS now points to ${BIGR_DEFAULT_TMP}"
  declare -x _JAVA_OPTIONS
  _JAVA_OPTIONS="-Djava.io.tmpdir='${BIGR_DEFAULT_TMP}'"
  export _JAVA_OPTIONS
elif [ "${_JAVA_OPTIONS:-}" == "/tmp" ]; then
  message WARNING "_JAVA_OPTIONS currently points to '/tmp'. This can lead to OS errors due to lack of space in /tmp"
  message WARNING "_JAVA_OPTIONS now points to ${BIGR_DEFAULT_TMP}"
  _JAVA_OPTIONS="-Djava.io.tmpdir='${BIGR_DEFAULT_TMP}'"
  export _JAVA_OPTIONS
else
  message INFO "_JAVA_OPTIONS -> ${_JAVA_OPTIONS}"
fi

if [ ! -f "${HOME}/.Renviron" ]; then
  message WARNING "${HOME}/.Renviron does not exists. This can lead to OS errors in R due to lack of space in /tmp"
  message WARNING "${HOME}/.Renviron was created with: TMP = '${BIGR_DEFAULT_TMP}'"
  echo -e "TMP = '${BIGR_DEFAULT_TMP}'" > "${HOME}/.Renviron"
fi

if [ ! -f "${HOME}/.condarc" ]; then
  message WARNING "${HOME}/.condarc does not exists. This can lead to OS errors in conda due to lack of space in /tmp"
  message WARNING "${HOME}/.condarc was created with: env_dir, pkgs_dir, and conda-build:root_dir overloaded"
  echo -e "envs_dir:\n\t- /mnt/beegfs/userdata/${USER}/anaconda/envs\npkgs_dir:\n\t- /mnt/beegfs/userdata/${USER}/anaconda/pkgs\nconda-build:\n\troot_dir: /mnt/beegfs/userdata/${USER}/conda-builds" > "${HOME}/.condarc"
fi


message INFO "sourced: environment.sh"
