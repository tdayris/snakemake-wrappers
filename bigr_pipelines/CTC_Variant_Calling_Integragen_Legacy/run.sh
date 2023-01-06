#!/usr/bin/env bash
set -e

PIPELINE_NAME=$(basename $(dirname "${0}"))

mkdir -vp data_{input,output}

bash "$(dirname ${0})/../common/bash/run.sh" --name ${PIPELINE_NAME} --use-singularity "$@"