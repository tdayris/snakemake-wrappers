#!/usr/bin/env bash

script_dir=$(dirname $0)
script_dir=$(readlink -e "${script_dir}")

if [[ ! "${PATH}" =~ conda ]]
then
    # If there is no available conda, use mine
    source "/mnt/beegfs/userdata/t_dayris/anaconda3/etc/profile.d/conda.sh"
fi

bash "${script_dir}/run.sh" -p demux output/multiqc.html "$@"
