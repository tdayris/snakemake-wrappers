#!/usr/bin/env bash

# This launcher is provided for piREST
# It call the usual run.sh file with parameters

genome=`cut -f4 design.tsv | grep -P "GRC" | sort | uniq | sed -P 's/"//g'`

if [ "${genome}" == "GRCh19" ]; then
  bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/cytoscan_eacon/run.sh ${genome}"
fi
