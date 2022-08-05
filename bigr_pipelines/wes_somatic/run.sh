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
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/wes_somatic"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH

SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config/config.hg38.yaml"
SNAKE_ARGS=("--cache estimate_igs_sureselect_v5")
PROFILE="slurm"
SUMMARY=""
GRAPH=""
ANMO="IGNORE"
CACHE_RULES=()

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph) GRAPH="${2}"; shift 2;;
    hg38|HG38) CONFIG_PATH="${PIPELINE_PATH}/config/config.hg38.yaml"; CACHE_RULES+=("msisensor_pro_scan_hg38");  CACHE_RULES+=("bwa_index_hg38"); shift;;
    mm10|MM10) CONFIG_PATH="${PIPELINE_PATH}/config/config.mm10.yaml"; CACHE_RULES+=("msisensor_pro_scan_mm10"); CACHE_RULES+=("bwa_index_mm10"); shift;;
    tmb|TMB) SNAKE_ARGS+=("--until tmb_only");;
    msi|MSI) SNAKE_ARGS+=("--until msi_only");;
    cnv|CNV|facets) SNAKE_ARGS+=("--until cnv_only");;
    map|mapping) SNAKE_ARGS+=("--until mapping_only");;
    fusions) SNAKE_ARGS+=("--until fusions_only");;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

if [ ! -f "config.yaml" ]; then
  message INFO "Missing config file."
  message INFO "Falling back to default arguments: ${CONFIG_PATH}"
  rsync -v "${CONFIG_PATH}" "config.yaml"
else
  message INFO "Config file already provided"
fi

if [ "${ANMO}" = "RUN" ]; then
  message INFO "ANMO-like filters are used on final TSV files"
  echo -e "\nANMO: true" >> "config.yaml"
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
