set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/clinics"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/cytoscan_eacon/Snakefile"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH
message INFO "Environment loaded"
SNAKE_ARGS=()

while [ "$#" -gt 0 ]; do
  case "${1}" in
    *) SNAKE_ARGS+=("${1}"); shift;;
  esac
done
message INFO "Environment loaded"

# Run pipeline
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"
snakemake -s "${PIPELINE_PATH}" --cache eacon_install eacon_databases --profile "${SNAKEMAKE_PROFILE_PATH}" "${SNAKE_ARGS[@]}" && message INFO "EaCoN successful" || error_handling "${LINENO}" 2 "Error while running Cytoscan pipeline"
