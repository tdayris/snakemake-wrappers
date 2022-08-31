#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX=$(readlink -e "$(dirname ${0})/../..")
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/messages.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/environment.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/CTC_Variant_Calling_Integragen_Legacy"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH

CWD=$(readlink -e ${PWD})
SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config/config.hg38.yaml"
SNAKE_ARGS=(
  "--use-envmodules" 
  "--use-singularity" 
  "--singularity-args '-B ${CWD}:/WORKDIR/${CWD} -B /mnt/beegfs/database/bioinfo/:/WORKDIR/mnt/beegfs/database/bioinfo/ -B ${CWD}/tmp/:/tmp/'"
)
PROFILE="slurm"
SUMMARY=""
GRAPH=""

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph|--dag) GRAPH="${2}"; shift 2;;
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

BASE_CMD="snakemake -s ${SNAKEFILE_PATH} --profile ${SNAKEMAKE_PROFILE_PATH}/${PROFILE} ${SNAKE_ARGS[*]}"
if [ "${SUMMARY}" != "" ]; then
  SUMMARY_CMD="${BASE_CMD} --summary > ${SUMMARY}"
  message CMD "${SUMMARY_CMD}"
  eval SUMMARY_CMD
elif [ "${GRAPH}" != "" ]; then
  RULEGRAPH_CMD="${BASE_CMD} --rulegraph | dot -Tpng > ${GRAPH}"
  message CMD "${RULEGRAPH_CMD}"
  eval ${SUMMARY_CMD}
else
  message CMD "${BASE_CMD}"
  eval ${BASE_CMD}
fi
