#!/usr/bin/env bash
set -euiop pipefail
shopt -s nullglob

OLD_VERSION="1.1.0"
NEW_VERSION="1.1.1"

find bigr_pipelines -type f -name "readme.md" | while read REAME; do
    sed -i "s|official-snakemake-wrappers/${OLD_VERSION}/bigr_pipelines|official-snakemake-wrappers/${NEW_VERSION}/bigr_pipelines|g" "${README}"
done