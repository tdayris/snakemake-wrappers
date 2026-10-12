#!/usr/bin/bash
set -euo 'pipefail'
shopt -s 'nullglob'

VERSION="0.5"

echo "Updating environment if needed"
echo "snakedeploy update-conda-envs '${PWD}/environment.yaml'"
echo "Building pin"
echo "snakedeploy pin-conda-envs '${PWD}/environment.yaml'"
echo "Building:"
echo "apptainer build '${PWD}/snakemake-conda-apptainer.sif' '${PWD}/snakemake-conda-apptainer.def'"
echo "Testing:"
echo "apptainer test '${PWD}/snakemake-conda-apptainer.sif'"
echo "Uploading:"
echo "apptainer push '${PWD}/snakemake-conda-apptainer.sif' 'oras://ghcr.io/tdayris/snakemake-conda-apptainer:${VERSION}'"
