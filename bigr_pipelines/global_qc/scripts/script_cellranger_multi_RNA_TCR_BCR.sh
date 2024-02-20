#!/bin/bash
set -e

#Script to launch cellranger multi analysis for single-cell RNA-seq data (+/-TCR +/-BCR)
#using: bash script_cellranger_multi_RNA_TCR_BCR.sh ${Thread} ${Memory} ${Sample_Name} ${CellRanger_param_file} ${pwd_folder}

#prepare environment
Thread=${1}
Memory=$((${2}/1000 - 1))
CR_sample=${3}
csv_config=${4}
pwd_folder=${5}

cd ${pwd_folder}
cd ./cellranger/

#if _loc file exist, CellRanger can't run (this file comes from a previous run):
[ -e ./${CR_sample}/_lock ] && rm ./${CR_sample}/_lock

#run cellranger
/tools/cellranger-7.2.0/bin/cellranger multi \
    --id=${CR_sample} \
    --csv=${csv_config} \
    --localcores=${Thread} \
    --localmem=${Memory}

#save results and log
mkdir -p ./web_summary/ ./csv_summary/
#cp ./${CR_sample}/outs/per_sample_outs/${CR_sample}/web_summary.html ./web_summary/${CR_sample}_web_summary_mqc.html
cp ./${CR_sample}/outs/per_sample_outs/${CR_sample}/metrics_summary.csv ./csv_summary/${CR_sample}_metrics_summary.csv
cp ./${CR_sample}/_log ../logs/cellranger/${CR_sample}_cellranger_multi.log

#clean folder
#rm -r ./${CR_sample}
