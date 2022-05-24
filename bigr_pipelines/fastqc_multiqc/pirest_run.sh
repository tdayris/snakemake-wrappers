#!/usr/bin/env bash

script_dir=$(dirname $0)
script_dir=$(realpath "${script_dir}")

bash "${script_dir}/run.sh" -p demux output/multiqc.html "$@"
