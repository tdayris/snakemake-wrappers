#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX=$(readlink -e "$(dirname ${0})/../..")
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/messages.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/environment.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/gene_set_enrichment"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH

SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"
SNAKE_ARGS=()
PROFILE="slurm"
SUMMARY=""
GRAPH=""

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph) GRAPH="${2}"; shift 2;;
    # hg19|HG19) CONFIG_PATH="${PIPELINE_PATH}/config.hg19.yaml"; shift;;
    hg38|HG38) CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"; shift;;
    # mm10|MM10) CONFIG_PATH="${PIPELINE_PATH}/config.mm10.yaml"; shift;;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

if [ ! -f "config.yaml" ]; then
  rsync -cv "${CONFIG_PATH}" "config.yaml"
else
  message INFO "Config file already provided"
fi

# Run pipeline
message CMD "conda_activate ${CONDA_ENV_PATH}"
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"

if [ "${SUMMARY}" != "" ]; then
  message CMD "snakemake -s ${SNAKEFILE_PATH} --configfile config.yaml --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} ${SNAKE_ARGS[*]} --summary > ${SUMMARY}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "config.yaml" --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" "${SNAKE_ARGS[@]}" --summary > "${SUMMARY}"
elif [ "${GRAPH}" != "" ]; then
  message CMD "snakemake -s ${SNAKEFILE_PATH} --configfile config.yaml --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} ${SNAKE_ARGS[*]} --rulegraph | dot -Tpng > ${GRAPH}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "config.yaml" --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" "${SNAKE_ARGS[@]}" --rulegraph | dot -Tpng > "${GRAPH}"
else
  message CMD "snakemake -s ${SNAKEFILE_PATH} --configfile config.yaml --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} ${SNAKE_ARGS[*]}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "config.yaml" --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" "${SNAKE_ARGS[@]}" && message INFO "GSEA successful" || error_handling "${LINENO}" 2 "Error while running GSEA pipeline"
fi