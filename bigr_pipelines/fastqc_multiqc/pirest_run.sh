#!/usr/bin/env bash

script_dir=$(dirname $0)
script_dir=$(readlink -e "${script_dir}")

params_demux=""

if [ -n "$(find input/*/archive/*/unaligned/Stats/ -name "Stats.json.zip" | head -1)" ]; then
    # Special case demux
    bash "${script_dir}/run.sh" -p demux output/multiqc.html "$@"
else
    echo "Stats not found"
    # Classic case, must fit demux-like behaviour
    bash "${script_dir}/run.sh" -p demux "$@"
    if [ -d "multiqc" ] ; then 
        mv --verbose multiqc output
    fi
fi