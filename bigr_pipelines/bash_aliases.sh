#!/bin/sh


# Activate main snakemake environment

function bigr_conda_activate () {
  source "$(conda info --base)/etc/profile.d/conda.sh" && 
  conda activate && 
  conda activate "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/env2/" ;
}

# Search a bigr project
function bigr_search_project () {
    ls -dlh /mnt/beegfs/scratch/bioinfo_core/*${1}*
    ls -dlh /mnt/isilon/data_bioinfo/Projets/*${1}* 
    ls -dlh /mnt/isilon/data_bioinfo/Projets/bioinfo_core/*${1}*
}

# Source this file to add new pipelines to your
# environment. No more loooong command lines,
# only the best.

alias bigr_fastq_multiqc="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/fastqc_multiqc/run.sh"
alias bigr_rnaseq="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh"
alias bigr_wes_somatic="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/wes_somatic/run.sh"
alias bigr_make_pon="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/make_pon/run.sh"
alias bigr_wes_germline="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_germline/run.sh"
alias bigr_snpeff_snpsift="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/snpeff_snpsift/run.sh"
alias bigr_ampliseq="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_ampliseq/run.sh"
alias bigr_sigprofiler_signatures="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/sigprofiler_signatures/run.sh"
alias bigr_oncoscan="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/oncoscan_eacon/run.sh"
alias bigr_cytoscan="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/cytoscan_eacon/run.sh"
alias bigr_ctc_legacy_integragen="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/CTC_Variant_Calling_Integragen_Legacy/run.sh"
alias bigr_bisulfite_ampliseq="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/bs_ampliseq/run.sh"
alias bigr_miraseq="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/miraseq/run.sh"
alias bigr_nfcore_chipseq="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/nfcore_chipseq/run.sh"


# Shortcuts for data managment
alias bigr_prepare_design="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/prepare_design.sh"
alias bigr_archive="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/archive.sh"

# Get a pipeline state and squeue data
alias bigr_status="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/status.sh"
alias bigr_my_jobs="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/my_jobs.sh"
alias bigr_all_jobs="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/all_jobs.sh"
alias bigr_post_pipeline="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/python/woops.py"

alias bigr_watch_my_jobs="watch --colors -n 10 bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/my_jobs.sh"
alias bigr_cancel_my_jobs="echo scancel $(squeue -u ${USER} | column -t |cut -f1 -d' ')"

# Multi process a bash command
alias bigr_multiprocess_bash="bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/bash/multiprocess_bash_command.sh"

# Goto functions
alias bigr_scratch_bioinfo="cd /mnt/beegfs/scratch/bioinfo_core/"
alias bigr_database="cd /mnt/beegfs/database/bioinfo/Index_DB/"
alias bigr_userdata="cd /mnt/beegfs/userdata/${USER}/"

alias bigr_meteo="curl https://fr.wttr.in/Villejuif"
alias bigr_best_employee="echo This user is the bast, by far: ${USER}"
alias bigr_list_users="who | cut -f 1 -d ' ' | sort | uniq -c"
