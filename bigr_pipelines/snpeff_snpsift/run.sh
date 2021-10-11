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
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/snpeff_snpsift"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH
message INFO "Environment loaded"

SNAKEFILE="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config.hg38.nochr.yaml"
SNAKE_ARGS=()
PROFILE="slurm"
SUMMARY=""
GRAPH=""

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph) GRAPH="${2}"; shift 2;;
    hg19chr|HG19chr) CONFIG_PATH="${PIPELINE_PATH}/config.hg19.yaml"; shift;;
    hg38chr|HG38chr) CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"; shift;;
    hg19|HG19) CONFIG_PATH="${PIPELINE_PATH}/config.hg19.nochr.yaml"; shift;;
    hg38|HG38) CONFIG_PATH="${PIPELINE_PATH}/config.hg38.nochr.yaml"; shift;;
    mm10|MM10) CONFIG_PATH="${PIPELINE_PATH}/config.mm10.yaml"; shift;;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

if [ ! -d "calls" ]; then
  error_handling "${LINENO}" 1 "VCF files must be in a directory called 'calls'"
fi

# Run pipeline
message CMD "conda_activate ${CONDA_ENV_PATH}"
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"

if [ "${SUMMARY}" != "" ]; then
  message CMD "snakemake -s ${SNAKEFILE} --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} --configfile ${CONFIG_PATH} ${SNAKE_ARGS[*]} --summary > ${SUMMARY}"
  snakemake -s "${SNAKEFILE}" --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" --configfile ${CONFIG_PATH} "${SNAKE_ARGS[@]}" --summary > "${SUMMARY}"
elif [ "${GRAPH}" != "" ]; then
  message CMD "snakemake -s ${SNAKEFILE} --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} --configfile ${CONFIG_PATH} ${SNAKE_ARGS[*]} --runegraph | dot -Tpng > ${GRAPH}"
  snakemake -s "${SNAKEFILE_PATH}" --configfile "config.yaml" --cache bwa_fixmate_meta_bwa_index bwa_index --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" "${SNAKE_ARGS[@]}" --rulegraph | dot -Tpng > "${GRAPH}"
else
  message CMD "snakemake -s ${SNAKEFILE} --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} --configfile ${CONFIG_PATH} ${SNAKE_ARGS[*]}"
  snakemake -s "${SNAKEFILE}" --profile "${SNAKEMAKE_PROFILE_PATH}/${PROFILE}" --configfile ${CONFIG_PATH} "${SNAKE_ARGS[@]}" && message INFO "SnpEff/SnpSift successful" || error_handling "${LINENO}" 2 "Error while running SnpEff/SnpSift pipeline"
fi
