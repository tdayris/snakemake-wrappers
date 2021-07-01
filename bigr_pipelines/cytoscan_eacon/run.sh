#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/messages.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/environment.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/cytoscan_eacon/Snakefile"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH
message INFO "Environment loaded"
SNAKE_ARGS=()
SUMMARY=""
PROFILE="slurm"
GRAPH=""

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph) GRAPH="${2}"; shift 2;;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

# Run pipeline
message CMD "conda_activate ${CONDA_ENV_PATH}"
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"

if [ "${SUMMARY}" != "" ]; then
  message CMD "snakemake -s ${PIPELINE_PATH} --cache eacon_install eacon_databases --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} ${SNAKE_ARGS[*]} --summary > ${SUMMARY}"
  snakemake -s "${PIPELINE_PATH}" --cache eacon_install eacon_databases --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" "${SNAKE_ARGS[@]}" --summary > "${SUMMARY}"
elif [ "${GRAPH}" != "" ]; then
  message CMD "snakemake -s ${PIPELINE_PATH} --cache eacon_install eacon_databases --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} ${SNAKE_ARGS[*]} --rulegraph | dot -Tpng > ${GRAPH}"
  snakemake -s "${PIPELINE_PATH}" --cache eacon_install eacon_databases --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" "${SNAKE_ARGS[@]}" --rulegraph | dot -Tpng > "${GRAPH}"
else
  message CMD "snakemake -s ${PIPELINE_PATH} --cache eacon_install eacon_databases --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} ${SNAKE_ARGS[*]}"
  snakemake -s "${PIPELINE_PATH}" --cache eacon_install eacon_databases --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" "${SNAKE_ARGS[@]}" && message INFO "EaCoN successful" || error_handling "${LINENO}" 2 "Error while running Cytoscan pipeline"
fi
