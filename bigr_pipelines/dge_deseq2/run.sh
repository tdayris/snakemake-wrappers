#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/messages.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/environment.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/dge_deseq2"
export PIPELINE_PATH

SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"
PROFILE="slurm"
SUMMARY=""
SNAKE_ARGS=()
GRAPH=""

while [ "$#" -gt 0 ]; do
  case "${1}" in
    hg19|HG19) CONFIG_PATH="${PIPELINE_PATH}/config.hg19.yaml"; shift;;
    hg38|HG38) CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"; shift;;
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph) GRAPH="${2}"; shift 2;;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

if [ ! -f "config.yaml" ]; then
  rsync -cv "${CONFIG_PATH}" "config.yaml"
else
  message INFO "Config file already provided"
fi

SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/${PROFILE}"

# Run pipeline
message CMD "conda_activate ${CONDA_ENV_PATH}"
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"

if [ "${SUMMARY}" != "" ]; then
  message INFO "snakemake -s ${SNAKEFILE_PATH} --configfile config.yaml --cache salmon_index --profile ${SNAKEMAKE_PROFILE_PATH} ${SNAKE_ARGS[*]} --summary > ${SUMMARY}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "config.yaml" --cache salmon_index --profile "${SNAKEMAKE_PROFILE_PATH}" "${SNAKE_ARGS[@]}" --summary > "${SUMMARY}"
elif [ "${GRAPH}" != "" ]; then
  message INFO "snakemake -s ${SNAKEFILE_PATH} --configfile config.yaml --cache salmon_index --profile ${SNAKEMAKE_PROFILE_PATH} ${SNAKE_ARGS[*]} --rulegraph | dot -Tpng > ${GRAPH}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "config.yaml" --cache salmon_index --profile "${SNAKEMAKE_PROFILE_PATH}" "${SNAKE_ARGS[@]}" --rulegraph | dot -Tpng > "${GRAPH}"
else
  message INFO "snakemake -s ${SNAKEFILE_PATH} --configfile config.yaml --cache salmon_index --profile ${SNAKEMAKE_PROFILE_PATH} ${SNAKE_ARGS[*]}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "config.yaml" --cache salmon_index --profile "${SNAKEMAKE_PROFILE_PATH}" "${SNAKE_ARGS[@]}"
fi
