set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

PROFILE="slurm";

while [ "$#" -gt 0 ]; do
  case "${1}" in
    -p|--profile) PROFILE="${2}"; shift 2;;
    -*|*) error_handling ${LINENO} 1 "Unknown arguments ${1}"; exit 1;;
  esac
done

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/${PROFILE}"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/fastqc_multiqc"
declare -x SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH SNAKEFILE_PATH
message INFO "Environment loaded"
message INFO "${SNAKEMAKE_PROFILE_PATH}"

# Run pipeline
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"
snakemake -s "${SNAKEFILE_PATH}" --profile "${SNAKEMAKE_PROFILE_PATH}"
