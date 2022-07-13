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
message INFO "Environment loaded"

README_PATH="${PIPELINE_PREFIX}/bigr_pipelines/README.md"
message INFO "Executing all tests written in ${README_PATH}, one after each other."

grep -Po "\-\>Test: \`\K[^\`]*?(?=\`)" "${README_PATH}" | while read COMMAND; do
  message CMD "${COMMAND}"
  ${COMMAND} > "test/all_tests.log" 2>&1
done
message INFO "Process over"