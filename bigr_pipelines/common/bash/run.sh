#!/usr/bin/env bash
set -euiop pipefail

# This adds several functions and variable in the environment
PIPELINE_PREFIX=$(readlink -e "$(dirname "${0}")/../../..")

COMMAND="source \"${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh\""
eval ${COMMAND}
message CMD "${COMMAND}"

COMMAND="source \"${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh\""
message CMD "${COMMAND}"
eval ${COMMAND}

# Define pipeline related variables
COMMAND="declare -x SNAKEMAKE_PROFILE_PATH=\"${PIPELINE_PREFIX}/bigr_pipelines/common/profiles\""
message CMD "${COMMAND}"
eval ${COMMAND}


COMMAND="declare -x WRAPPERS_PATH=\"${PIPELINE_PREFIX}\""
message CMD "${COMMAND}"
eval ${COMMAND}

COMMAND="export SNAKEMAKE_PROFILE_PATH WRAPPERS_PATH"
message CMD "${COMMAND}"
eval ${COMMAND}


# Default IO directories
CWD=$(readlink -e "${PWD}")
iodirectories "${CWD}" "data_input"
iodirectories "${CWD}" "data_output"

# Flamingo databases path
FLAMINGO_DB="/mnt/beegfs/database/bioinfo/"

# Default snakemake arguments
SNAKE_ARGS=(
  "--wrapper-prefix" "${WRAPPERS_PATH}/"
  "--singularity-args '-B ${SHARED_CONDA_INDSTALL}:${SHARED_CONDA_INDSTALL} -B ${CONDA_ENV_PATH}:${CONDA_ENV_PATH} -B ${SNAKEMAKE_OUTPUT_CACHE}:${SNAKEMAKE_OUTPUT_CACHE} -B ${CONDA_CACHE_PATH}:${CONDA_CACHE_PATH} -B ${CWD}:${CWD} -B ${FLAMINGO_DB}:${FLAMINGO_DB} -B ${CWD}/tmp/:/tmp/'"
)
STEPS=()
PROFILE="slurm"
SUMMARY=""
GRAPH=""
NAME=""
CONFIG_PATH=""

# Command line parser
while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    --summary) SUMMARY="${2}"; shift 2;;
    --rulegraph|--dag) GRAPH="${2}"; shift 2;;
    hg19|HG19|GRCh37) CONFIG_PATH="${PIPELINE_PATH}/config/config.hg19.yaml"; message WARNING "Some operations are not available for hg19 genome"; shift;;
    hg38|HG38|GRCh38) CONFIG_PATH="${PIPELINE_PATH}/config/config.hg38.yaml"; shift;;
    mm10|MM10|GRCm38) CONFIG_PATH="${PIPELINE_PATH}/config/config.mm10.yaml"; message WARNING "Some operations are not available for mice datasets"; shift;;
    DESeq2|deseq2|DGE|dge) STEPS+=("dge"); message INFO "DGE is in the expected result list"; shift;;
    salmon|Salmon|quant) STEPS+=("quant"); message INFO "Quantification is in the expected result list"; shift;;
    fusions|fusion) STEPS+=("fusions"); message INFO "Fusions are in the expected result list"; shift;;
    qc|QC) STEPS+=("qc"); message INFO "Quality Controls are in the expected result list"; shift;;
    immu|deconv) STEPS+=("immunedeconv"); message INFO "Immune Deconvolution is in the expected result list"; shift;;
    gsea|clusterprofiler) STEPS+=("gsea"); message INFO "GSEA are in the expected result list"; shift;;
    tmb|TMB) STEPS+=("tmb"); message INFO "TMB are in the expected result list"; shift;;
    msi|MSI) STEPS+=("msi"); message INFO "MSI are in the expected result list"; shift;;
    cnv|CNV|facets) STEPS+=("cnv"); message INFO "CNV are in the expected result list"; shift;;
    map|mapping) STEPS+=("aln"); message INFO "Alignments are in the expected result list"; shift;;
    vc|variant|variants|calling) STEPS+=("variant_calling"); message INFO "Variant callings are in the expected result list"; shift;;
    --name) NAME="${2}"; shift 2;;
    -h|--help) message INFO "Please see official documentation at: https://github.com/tdayris/snakemake-wrappers"; exit 0; shift;;
    --snake-help) SNAKE_ARGS+=("--help"); shift;;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Command line parsed"


message INFO "Working with pipeline ${NAME}"
COMMAND="declare -x PIPELINE_PATH=\"${PIPELINE_PREFIX}/bigr_pipelines/${NAME}\""
message CMD "${COMMAND}"
eval ${COMMAND}

COMMAND="export PIPELINE_PATH"
message CMD "${COMMAND}"
eval ${COMMAND}


# Default snakefile path
if [ -f "${PIPELINE_PATH}/Snakefile" ]; then
  SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
  message INFO "Snakefile found: ${SNAKEFILE_PATH}"
else
  message ERROR "Could not find Snakefile in: ${PIPELINE_PATH}"
  exit 1
fi

# Default config file path

if [ ! -f "${CONFIG_PATH}" ]; then
  if [ -f "${PIPELINE_PATH}/config/config.yaml" ]; then
    CONFIG_PATH="${PIPELINE_PATH}/config/config.yaml"
  elif [ -f "${PIPELINE_PATH}/config/config.hg38.yaml" ]; then
    CONFIG_PATH="${PIPELINE_PATH}/config/config.hg38.yaml"
  elif [ -f "${PIPELINE_PATH}/config.hg38.yaml" ]; then
    CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"
  fi
fi


# If config is not available in local repository, then add it !
if [ ! -f "config.yaml" ]; then
  if [ ! -f "${CONFIG_PATH}" ]; then
    message ERROR "Config file does not exist at ${CONFIG_PATH}. The ${NAME} pipeline may not be available for this genome."
  fi
  message INFO "Config file not found, falling back to default arguments."
  COMMAND="rsync --verbose ${CONFIG_PATH} config.yaml"
  message CMD "${COMMAND}"
  eval ${COMMAND}
else
  message INFO "Config file already provided"
fi

# If setps are defined in command line, the activate them
if [ ${#STEPS[@]} -gt 0 ]; then
  message INFO "Activating expected steps"
  for STEP in "${STEPS[@]}"; do
    COMMAND="sed -i 's/  ${STEP}: false/  ${STEP}: true/g' config.yaml"
    message CMD "${COMMAND}"
    eval ${COMMAND}
  done
else
  message INFO "No step to activate. Running default pipeline behaviour"
fi

# Run pipeline
# Activate singularity
CMD="module load singularity/3.6.3"
message CMD "${CMD}"
eval ${CMD} && message INFO "Singularity loaded in case your pipeline requires it."

# Activate conda
message CMD "conda_activate ${CONDA_ENV_PATH}"
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"

# Launch snakemake
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

message INFO "Gathering cluster usage information"
conda_activate "/mnt/beegfs/pipelines/unofficial-snakemake-wrappers/shared_install/Woops" || message ERROR "Could not read cluster usage statistics, but pipeline was OK and results can be trusted."
CMD="python3 ${PIPELINE_PREFIX}/bigr_pipelines/common/python/woops.py > logs/woops.log 2>&1"
message CMD "${CMD}"
eval ${CMD} || message ERROR "Could not read cluster usage statistics, but pipeline was OK and results can be trusted."
message INFO "See run info in logs."

message INFO "Process over."