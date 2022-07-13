#!/usr/bin/env bash

script_dir=$(dirname $0)
script_dir=$(readlink -e "${script_dir}")

bash "${script_dir}/run.sh" -p demux output/multiqc.html "$@"
