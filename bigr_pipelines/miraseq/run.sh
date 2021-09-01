#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/messages.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/environment.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/miraseq"

PROFILE="slurm";
CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml";
SUMMARY="";
SNAKE_ARGS=()
GRAPH="";

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

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/${PROFILE}"
declare -x SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH SNAKEFILE_PATH
message INFO "Environment loaded"
message INFO "${SNAKEMAKE_PROFILE_PATH}"

if [ ! -f "config.yaml" ]; then
  rsync -cv "${CONFIG_PATH}" "config.yaml"
else
  message INFO "Config file already provided"
fi

# Run pipeline
message CMD "conda_activate ${CONDA_ENV_PATH}"
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"

if [ "${SUMMARY}" != "" ]; then
  message CMD "snakemake -s ${SNAKEFILE_PATH} --configfile ${CONFIG_PATH} --profile ${SNAKEMAKE_PROFILE_PATH} ${SNAKE_ARGS[*]} --summary > ${SUMMARY}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "${CONFIG_PATH}" --profile "${SNAKEMAKE_PROFILE_PATH}" "${SNAKE_ARGS[@]}" --summary > "${SUMMARY}"
elif [ "${GRAPH}" != "" ]; then
  message CMD "snakemake -s ${SNAKEFILE_PATH} --configfile ${CONFIG_PATH} --profile ${SNAKEMAKE_PROFILE_PATH} ${SNAKE_ARGS[*]} --rulegraph | dot -Tpng > ${GRAPH}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "${CONFIG_PATH}" --profile "${SNAKEMAKE_PROFILE_PATH}" "${SNAKE_ARGS[@]}" --rulegraph | dot -Tpng > "${GRAPH}"
else
  message CMD "snakemake -s ${SNAKEFILE_PATH} --configfile ${CONFIG_PATH} --profile ${SNAKEMAKE_PROFILE_PATH} ${SNAKE_ARGS[*]}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "${CONFIG_PATH}" --profile "${SNAKEMAKE_PROFILE_PATH}" "${SNAKE_ARGS[@]}" && message INFO "MIRA-Seq pipeline successful" || error_handling "${LINENO}" 2 "Error while running MIRA-Seq pipeline"
fi
