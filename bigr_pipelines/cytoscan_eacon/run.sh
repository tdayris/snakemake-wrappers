set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/functions.sh"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/profile/clinics"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/meta/bio/eacon_cytoscan"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH

# Run pipeline
conda activate "${CONDA_ENV_PATH}" || error_handling "${LINENO}" 1 "Could not activate conda environment"
snakemake -s "${PIPELINE_PATH}/Snakefile" --profile "${SNAKEMAKE_PROFILE_PATH}" || error_handling "${LINENO}" 2 "Error while running Cytoscan pipeline"
