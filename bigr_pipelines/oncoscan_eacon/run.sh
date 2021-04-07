set -e

declare -x CONDA_ENV_PATH="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/env/"
declare -x SNAKEMAKE_PROFILE_PATH="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/profile/clinics"
declare -x SNAKEMAKE_OUTPUT_CACHE="/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/cache"
declare -x PIPELINE_PATH="/mnt/beegfs/pipelines/snakemake-wrappers/meta/bio/eacon_oncoscan"
export CONDA_ENV_PATH SNAKEMAKE_PROFILE_PATH SNAKEMAKE_OUTPUT_CACHE

conda activate "${CONDA_ENV_PATH}"
snakemake -s "${PIPELINE_PATH}/Snakefile" --profile "${SNAKEMAKE_PROFILE_PATH}"
