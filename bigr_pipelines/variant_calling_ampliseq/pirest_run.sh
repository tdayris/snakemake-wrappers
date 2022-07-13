#!/usr/bin/env bash

# This launcher is provided for piREST
# It call the usual run.sh file with parameters

genome=`cut -f4 design.tsv | grep -P "GRC" | sort | uniq | sed -P 's/"//g'`

bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_ampliseq/run.sh ${genome}"
