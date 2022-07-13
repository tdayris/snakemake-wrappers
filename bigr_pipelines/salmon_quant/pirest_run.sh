#!/usr/bin/env bash

# This launcher is provided for piREST
# It call the usual run.sh file with parameters

genome=`cut -f4 design.tsv | grep -P "GRC" | sort | uniq | sed 's/"//g'`
script_dir=$(dirname $0)
script_dir=$(readlink -e "${script_dir}")

bash "${script_dir}/run.sh" "${genome}" --nt --rerun-incomplete
