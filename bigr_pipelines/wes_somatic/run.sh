#!/usr/bin/env bash
set -e

PIPELINE_NAME=$(basename $(dirname "${0}"))
ADDITIONAL_OPTS="--cache estimate_igs_sureselect_v5"

bash "$(dirname ${0})/../common/bash/run.sh" --name ${PIPELINE_NAME} "$@"