#!/bin/bash

set -eu

# Authors: Thibault Dayris
# Copyright: Copyright 2021, Thibault Dayris
# Email: thibault.dayris@gustaveroussy.fr
# License: MIT

# This adds several functions and variable in the environment
PIPELINE_PREFIX="/mnt/beegfs/pipelines/snakemake-wrappers"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/messages.sh"
source "${PIPELINE_PREFIX}/bigr_pipelines/common/bash/environment.sh"

declare -x CONDA_ENV_FILE="${PIPELINE_PREFIX}/bigr_pipelines/common/snakemake.yaml"
declare -x ENV_DESTINATION="${PIPELINE_PREFIX}/bigr_pipelines/common/env/"
export CONDA_ENV_FILE ENV_DESTINATION
message INFO "Environment updated"

message WARNING "You are about to change BiGR whole conda environment. This could impact other users."
message WARNING "Are you sure you want to do this?"
message WARNING "Did you ask Marc or Thibault?"
select yn in "Yes" "No"; do 
  case $REPLY in
    Yes ) mamba env create --force --file "${CONDA_ENV_FILE}" --prefix "${ENV_DESTINATION}" && message INFO "Installation successful" || error_handling "${LINENO}" 2 "Installation failed :-("; break;;
    No ) message INFO "No installation/update was performed."; break;;
    * ) message ERROR "Please answer Yes or No";;
  esac
done