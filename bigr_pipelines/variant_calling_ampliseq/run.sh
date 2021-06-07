set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

# Define pipeline related variables
declare -x SNAKEMAKE_PROFILE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/common/profiles/slurm"
declare -x PIPELINE_PATH="${PIPELINE_PREFIX}/bigr_pipelines/variant_calling_ampliseq"
export SNAKEMAKE_PROFILE_PATH PIPELINE_PATH

SNAKEFILE_PATH="${PIPELINE_PATH}/Snakefile"
CONFIG_PATH="${PIPELINE_PATH}/config.hg19.yaml"
if [ "${1}" = "hg38" ]; then
  CONFIG_PATH="${PIPELINE_PATH}/config.hg38.yaml"
  message WARNING "HG38 has not been tested."
fi
message INFO "Environment loaded"

if [ -f "config.yaml" ]; then
  if [ ! -f "config_variant_calling_ampliseq.yaml" ]; then
    cp -v "config.yaml" "config_variant_calling_ampliseq.yaml"
  fi
fi

if [ ! -f "config_variant_calling_ampliseq.yaml" ]; then
  message INFO "Missing config file, falling back to default arguments"
  cp "${CONFIG_PATH}" "config_variant_calling_ampliseq.yaml"
else
  message INFO "Config file already provided"
fi

# Run pipeline
conda_activate "${CONDA_ENV_PATH}" && message INFO "Conda loaded" || error_handling "${LINENO}" 1 "Could not activate conda environment"
snakemake -s "${SNAKEFILE_PATH}" --configfile "config_variant_calling_ampliseq.yaml" --cache bwa_fixmate_meta_bwa_index bwa_index --profile "${SNAKEMAKE_PROFILE_PATH}" && message INFO "Variant calling successful" || error_handling "${LINENO}" 2 "Error while running variant calling pipeline"
