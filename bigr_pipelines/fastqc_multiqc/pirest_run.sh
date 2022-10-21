#!/usr/bin/env bash

script_dir=$(dirname $0)
script_dir=$(readlink -e "${script_dir}")

if [[ ! "${PATH}" =~ conda ]]
then
    # If there is no available conda, use mine
    source "/mnt/beegfs/userdata/t_dayris/anaconda3/etc/profile.d/conda.sh"
fi

params_demux="-p demux"

if [ -n "$(find input/*/archive/*/unaligned/Stats/ -name "Stats.json.zip" | head -1)" ]; then
    # Special case demux
    bash "${script_dir}/run.sh" output/multiqc.html ${params_demux} "$@"
else
    echo "Stats not found"
    # Classic case, must fit demux-like behaviour
    bash "${script_dir}/run.sh" ${params_demux} "$@"
    if [ -d "multiqc" ] ; then 
        mv --verbose --force --update multiqc/* output
    fi
fi
