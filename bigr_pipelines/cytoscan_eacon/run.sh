set -e

# This adds several functions and variable in the environment
source "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/base_runner.sh"

declare -x SNAKEMAKE_PROFILE_PATH="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/profile/clinics"
declare -x PIPELINE_PATH="/mnt/beegfs/pipelines/snakemake-wrappers/meta/bio/eacon_cytoscan"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH

conda activate "${CONDA_ENV_PATH}" || error_handling "${LINENO}" 1 "Could not activate conda environment"
snakemake -s "${PIPELINE_PATH}/Snakefile" --profile "${SNAKEMAKE_PROFILE_PATH}" || error_handling "${LINENO}" 2 "Error while running Cytoscan pipeline"
