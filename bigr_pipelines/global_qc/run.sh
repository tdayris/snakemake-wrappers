#!/usr/bin/env bash
set -e

PIPELINE_NAME=$(basename $(dirname "${0}"))

bash "$(dirname ${0})/../common/bash/run.sh" --name ${PIPELINE_NAME} "$@"