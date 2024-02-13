#!/bin/bash
set -e

#Script to launch cellranger atac analysis for single-cell ATAC-seq data
#using: bash script_cellranger_atac_ATAC.sh ${Thread} ${Memory} ${Sample_Name} ${Library_Name} ${reference} ${path_to_fastq_files} ${pwd_folder}

#prepare environment
Thread=${1}
Memory=$((${2}/1000 - 1))
CR_sample=${3}
library_names=${4}
reference=${5}
path_fastqs=${6}
pwd_folder=${7}

cd ${pwd_folder}
cd ./cellranger/

#if _loc file exist, CellRanger can't run (this file comes from a previous run):
[ -e ./${CR_sample}/_lock ] && rm ./${CR_sample}/_lock

#run cellranger
/tools/cellranger-atac-2.1.0/cellranger-atac count \
    --id=${CR_sample} \
    --sample=${library_names} \
    --reference=${reference} \
    --fastqs=${path_fastqs} \
    --localcores=${Thread} \
    --localmem=${Memory}
#reference=/mnt/beegfs/database/bioinfo/cellranger-atac/2.1.0/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 # todo !!!!

#save results and log
mkdir -p ./web_summary/ ./csv_summary/
#cp ./${CR_sample}/outs/web_summary.html ./web_summary/${CR_sample}_web_summary_mqc.html
cp ./${CR_sample}/outs/summary.csv ./csv_summary/${CR_sample}_metrics_summary.csv
cp ./${CR_sample}/_log ../logs/cellranger/${CR_sample}_cellranger_atac.log

#clean folder
#rm -r ./${CR_sample}
