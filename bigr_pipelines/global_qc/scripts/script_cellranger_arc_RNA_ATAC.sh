#!/bin/bash
set -e

#Script to launch cellranger arc analysis for single-cell ATAC-seq + RNA-seq data
#using: bash script_cellranger_arc_RNA_ATAC.sh ${Thread} ${Memory} ${Sample_Name} ${CellRanger_param_file} ${pwd_folder}

#prepare environment
Thread=${1}
Memory=$((${2}/1000 - 1))
CR_sample=${3}
csv_config=${4}
reference=${5}
pwd_folder=${6}

cd ${pwd_folder}
cd ./cellranger/

#if _loc file exist, CellRanger can't run (this file comes from a previous run):
[ -e ./${CR_sample}/_lock ] && rm ./${CR_sample}/_lock

#run cellranger
/tools/cellranger-arc-2.0.2/bin/cellranger-arc count \
    --id=${CR_sample} \
    --reference=${reference} \
    --libraries=${csv_config} \
    --localcores=${Thread} \
    --localmem=${Memory}

#save results and log
mkdir -p ./web_summary/ ./csv_summary/
#cp ./${CR_sample}/outs/web_summary.html ./web_summary/${CR_sample}_web_summary_mqc.html
cp ./${CR_sample}/outs/summary.csv ./csv_summary/${CR_sample}_metrics_summary.csv
cp ./${CR_sample}/_log ../logs/cellranger/${CR_sample}_cellranger_arc.log

#clean folder
#rm -r ./${CR_sample}


