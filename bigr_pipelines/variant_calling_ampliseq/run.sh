set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/slurm"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/variant_calling_ampliseq"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH
message INFO "Environment loaded"

SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
if [ "${1}" = "hg19" ]; then
  SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile_hg19.smk"
fi

# Run pipeline
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"
snakemake -s "${SNAKEFILE_PATH}" --cache eacon_install eacon_databases --profile "${SNAKEMAKE_PROFILE_PATH}" && message INFO "Variant calling successful" || error_handling "${LINENO}" 2 "Error while running variant calling pipeline"
