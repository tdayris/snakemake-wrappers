#!/usr/bin/bash

set -euo 'pipefail'
shopt -s 'nullglob'

PIPELINE_DIR=$(dirname $(readlink -e "${0}"))
MAIN_RUN_SCRIPT=$(echo -e "$(dirname "${PIPELINE_DIR}")/run.py")

eval "$(conda shell hook --shell bash)"
conda activate snakemake2

python "${MAIN_RUN_SCRIPT}" --pipeline-name "genomes_download" --profile "local"
