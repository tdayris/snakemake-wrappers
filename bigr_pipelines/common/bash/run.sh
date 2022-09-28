#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX=$(readlink -e "$(dirname ${0})/../../..")
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/messages.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
# shellcheck source=/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/environment.sh
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/${1}"
declare -x WRAPPERS_PATH=$(readlink -e "${PIPELINE_PATH}/../../")
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH WRAPPERS_PATH

if [ -f "${PIPELINE_PATH}/Snakefile" ]; then
    SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
else
    message ERROR "Could not find Snakefile in: ${PIPELINE_PATH}"
    exit 1
fi

if [ -f "${PIPELINE_PATH}/config/config.hg38.yaml" ]; then
    CONFIG_PATH="${PIPELINE_PATH}/config/config.hg38.yaml"
elif [ -f "${PIPELINE_PATH}/config.hg38.yaml" ]; then
    CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"
fi

SNAKE_ARGS=("--wrapper-prefix" "${WRAPPERS_PATH}/")
STEPS=()
PROFILE="slurm"
SUMMARY=""
GRAPH=""

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph|--dag) GRAPH="${2}"; shift 2;;
    hg19|HG19|GRCh37) CONFIG_PATH="${PIPELINE_PATH}/config/config.hg19.yaml"; message WARNING "Some operations are not available for hg19 genome"; shift;;
    hg38|HG38|GRCh38) shift;;
    mm10|MM10|GRCm38) CONFIG_PATH="${PIPELINE_PATH}/config/config.mm10.yaml"; message WARNING "Some operations are not available for mice datasets"; shift;;
    DESeq2|deseq2|DGE|dge) STEPS+=("dge"); message INFO "DGE is in the expected result list"; shift;;
    salmon|Salmon|quant) STEPS+=("quant"); message INFO "Quantification is in the expected result list"; shift;;
    fusions|fusion) STEPS+=("fusions"); message INFO "Fusions are in the expected result list"; shift;;
    qc|QC) STEPS+=("qc"); message INFO "Quality Controls is in the expected result list"; shift;;
    immu|deconv) STEPS+=("immunedeconv"); message INFO "Immune Deconvolution is in the expected result list"; shift;;
    gsea|clusterprofiler) STEPS+=("gsea"); message INFO "GSEA is in the expected result list"; shift;;
    rnaseq) shift;;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

if [ ! -f "config.yaml" ]; then
  if [ ! -f "${CONFIG_PATH}" ]; then
    message ERROR "Config file does not exist at ${CONFIG_PATH}. The ${1} pipeline may not be available for this genome."
  fi
  COMMAND="rsync --verbose ${CONFIG_PATH} config.yaml"
  message CMD "${COMMAND}"
  eval ${COMMAND}
else
  message INFO "Config file already provided"
fi

message INFO "Activating expected steps if available in the pipeline"
for STEP in "${STEPS[@]}"; do
  COMMAND="sed -i 's/^  ${STEP}: false$/  ${STEP}: true/g' config.yaml"
  message CMD "${COMMAND}"
  eval ${COMMAND}
done

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
