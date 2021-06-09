set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/slurm"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/salmon_quant"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH

SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"
if [ "${1}" = "hg19" ]; then
  CONFIG_PATH="${PIPELINE_PATH}/config.hg19.yaml"
  message INFO "Using HG19"
fi
message INFO "Environment loaded"

if [ ! -f "config.yaml" ]; then
  rsync -cv "${CONFIG_PATH}" "config.yaml"
else
  message INFO "Config file already provided"
fi

# Run pipeline
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"
snakemake -s "${SNAKEFILE_PATH}" --configfile "config.yaml" --cache salmon_index --profile "${SNAKEMAKE_PROFILE_PATH}" && message INFO "Variant calling successful" || error_handling "${LINENO}" 2 "Error while running Salmon quant pipeline"
