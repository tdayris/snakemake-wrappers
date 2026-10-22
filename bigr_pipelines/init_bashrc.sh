#!/bin/bash
set -euo 'pipefail'
shopt -s 'nullglob'

. /etc/bash.bashrc
bash /etc/bash.bashrc
eval "$(/root/miniforge3/bin/conda 'shell.bash' 'hook')" 
conda info --envs
conda activate /root/miniforge3/envs/snakemake
exec "$@"
