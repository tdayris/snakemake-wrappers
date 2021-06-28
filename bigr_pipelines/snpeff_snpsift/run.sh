set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/slurm"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/snpeff_snpsift"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH
message INFO "Environment loaded"

SNAKEFILE="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"
SNAKE_ARGS=()

while [ "$#" -gt 0 ]; do
  case "${1}" in
    hg19|HG19) CONFIG_PATH="${PIPELINE_PATH}/config.hg19.yaml"; shift;;
    hg38|HG38) CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"; shift;;
    hg19nochr) CONFIG_PATH="${PIPELINE_PATH}/config.hg19.nochr.yaml"; shift;;
    hg38nochr) CONFIG_PATH="${PIPELINE_PATH}/config.hg38.nochr.yaml"; shift;;
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

if [ ! -d "calls" ]; then
  error_handling "${LINENO}" 1 "VCF files must be in a directory called 'calls'"
fi

# Run pipeline
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"
snakemake -s "${SNAKEFILE}" --profile "${SNAKEMAKE_PROFILE_PATH}" --configfile ${CONFIG_PATH} "${SNAKE_ARGS[@]}" && message INFO "SnpEff/SnpSift successful" || error_handling "${LINENO}" 2 "Error while running SnpEff/SnpSift pipeline"
