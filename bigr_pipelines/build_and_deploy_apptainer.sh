#!/usr/bin/bash
set -euo 'pipefail'
shopt -s 'nullglob'

VERSION="0.14"

echo -e "\CONDA\n"
echo "Updating environment if needed"
echo "snakedeploy update-conda-envs '${PWD}/environment.yaml'"
echo "Building pin"
echo "snakedeploy pin-conda-envs '${PWD}/environment.yaml'"

echo -e "\nAPPTAINER\n"
echo "Building:"
echo "apptainer build '${PWD}/snakemake-conda-apptainer.sif' '${PWD}/snakemake-conda-apptainer.def'"
echo "Testing:"
echo "apptainer test '${PWD}/snakemake-conda-apptainer.sif'"
echo "Uploading:"
echo "apptainer push '${PWD}/snakemake-conda-apptainer.sif' 'oras://ghcr.io/tdayris/snakemake-conda-apptainer:${VERSION}'"

echo -e "\nDOCKER\n"
echo "Login if needed:"
echo "docker login ghcr.io -u tdayris"
echo "Translating apptainer to Dockerfile"
echo "spython recipe --parser singularity --writer docker snakemake-conda-apptainer-python.def | sed '$ d' > Dockerfile"
echo "Build docker image from dockerfile"
echo "docker build -t snakemake-conda-apptainer:${VERSION} -t snakemake-conda-apptainer:latest ."
echo "Define release tag"
echo "docker tag snakemake-conda-apptainer:${VERSION} ghcr.io/tdayris/snakemake-conda-apptainer:${VERSION}"
echo "Push to github"
echo "docker push ghcr.io/tdayris/snakemake-conda-apptainer:${VERSION}"
