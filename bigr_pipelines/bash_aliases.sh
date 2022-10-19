#!/bin/sh
SNAKEMAKE_WRAPPER_VERSION="1.1.0"
SNAKEMAKE_WRAPPER_PREFIX=$(readlink -e "${0}/../../../")

# Activate main snakemake environment

function bigr_conda_activate () {
  source "$(conda info --base)/etc/profile.d/conda.sh" && 
  conda activate && 
  conda activate "${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/bigr_snakemake/" ;
}

# Search a bigr project
function bigr_dm_search_project () {
    ls -dlh /mnt/beegfs/scratch/bioinfo_core/*${1}*
    ls -dlh /mnt/isilon/data_bioinfo/Projets/*${1}* 
    ls -dlh /mnt/isilon/data_bioinfo/Projets/bioinfo_core/*${1}*
}

# Deploy new version on Flamingo
function bigr_dm_deploy_new_version () {
  NEW_VERSION="${1}"
  source "${SNAKEMAKE_WRAPPER_PREFIX}/bigr_pipeline/common/bash/messages.sh"

  if [ -d "${SNAKEMAKE_WRAPPER_PREFIX}/${NEW_VERSION}" ]; then
    message ERROR "The suggested version already exists, please provide a new number or update the available one."
  else
    message INFO "Building new version directory"
    COMMAND="mkdir --parents --verbose \"${SNAKEMAKE_WRAPPER_PREFIX}/${NEW_VERSION}\""
    message CMD "${COMMAND}"
    eval ${COMMAND}

    message INFO "Clonig git repository"
    COMMAND="git clone https://github.com/tdayris/snakemake-wrappers.git \"${SNAKEMAKE_WRAPPER_PREFIX}/${NEW_VERSION}\""
    message CMD "${COMMAND}"
    eval ${COMMAND}

    message INFO "Reaching last commit on Unofficial"
    COMMAND="git checkout Unofficial"
    message CMD "${COMMAND}"
    eval ${COMMAND}

    message INFO "Sourcing new version for today"
    COMMAND="source \"${SNAKEMAKE_WRAPPER_PREFIX}/${NEW_VERSION}/bigr_pipeline/bash_aliases.sh\""
    message CMD "${COMMAND}"
    eval ${COMMAND}

    message WARNING "Do not forget to add the above command to your bashrc, and remove the old one!"
    message INFO "Process over"
  fi
}

# Source this file to add new pipelines to your
# environment. No more loooong command lines,
# only the best.

alias bigr_pipeline_fastq_multiqc="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/fastqc_multiqc/run.sh"
alias bigr_pipeline_rnaseq="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/rnaseq/run.sh"
alias bigr_pipeline_wes_somatic="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/wes_somatic/run.sh"
alias bigr_pipeline_make_pon="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/make_pon/run.sh"
alias bigr_pipeline_wes_germline="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/variant_calling_germline/run.sh"
alias bigr_pipeline_snpeff_snpsift="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/snpeff_snpsift/run.sh"
alias bigr_pipeline_ampliseq="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/variant_calling_ampliseq/run.sh"
alias bigr_pipeline_sigprofiler_signatures="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/sigprofiler_signatures/run.sh"
alias bigr_pipeline_oncoscan="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/oncoscan_eacon/run.sh"
alias bigr_pipeline_cytoscan="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/cytoscan_eacon/run.sh"
alias bigr_pipeline_ctc_legacy_integragen="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/CTC_Variant_Calling_Integragen_Legacy/run.sh"
alias bigr_pipeline_bisulfite_ampliseq="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/bs_ampliseq/run.sh"
alias bigr_pipeline_miraseq="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/miraseq/run.sh"
alias bigr_pipeline_nfcore_chipseq="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/nfcore_chipseq/run.sh"


# Shortcuts for data managment
alias bigr_dm_prepare_design="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/common/bash/prepare_design.sh"
alias bigr_dm_archive="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/common/bash/archive.sh"

# Get a pipeline state and squeue data
alias bigr_pipeline_status="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/common/bash/status.sh"
alias bigr_slurm_my_jobs="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/common/bash/my_jobs.sh"
alias bigr_slurm_all_jobs="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/common/bash/all_jobs.sh"
alias bigr_pipeline_post_process="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/common/python/woops.py"

alias bigr_slurm_watch_my_jobs="watch --colors -n 10 bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/common/bash/my_jobs.sh"
alias bigr_slurm_cancel_my_jobs="echo scancel $(squeue -u ${USER} | column -t |cut -f1 -d' ')"
alias bigr_utils_list_users="who | cut -f 1 -d ' ' | sort | uniq -c"

# Multi process a bash command
alias bigr_utils_multiprocess_bash="bash ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/common/bash/multiprocess_bash_command.sh"
alias bigr_utils_update_dev_git="cd /mnt/beegfs/userdata/${USER}/snakemake-wrappers/ && git pull && cd -"
alias bigr_utils_update_prod_git="cd ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION} && git pull && cd -"

# Goto functions
alias bigr_cd_scratch_bioinfo="cd /mnt/beegfs/scratch/bioinfo_core/"
alias bigr_cd_database="cd /mnt/beegfs/database/bioinfo/Index_DB/"
alias bigr_cd_userdata="cd /mnt/beegfs/userdata/${USER}/"
alias bigr_cd_pipeline="cd ${SNAKEMAKE_WRAPPER_PREFIX}/${SNAKEMAKE_WRAPPER_VERSION}/"

# Today: Fun with flags
alias bigr_utils_meteo="curl https://fr.wttr.in/Villejuif"
alias bigr_utils_best_employee="echo This user is the best, by far: ${USER}"
alias bigr_utils_version="echo Using Unofficial Snakemake-wrapper, version: ${SNAKEMAKE_WRAPPER_VERSION}"
