#!/bin/bash
set -euo 'pipefail'
shopt -s 'nullglob'

export SNAKEMAKE_GOOGLEBATCH_PROJECT="canceropole-infra"
export GOOGLE_CLOUD_PROJECT="canceropole-infra"
export SNAKEMAKE_GOOGLEBATCH_REGION="europe-west9"

IN_CREDENTIALS="${HOME}/.config/gcloud/application_default_credentials.json"

docker run -it \
    -v "${IN_CREDENTIALS}":/config/gcloud/application_default_credentials.json:ro \
    -v .:/app \
    -e GOOGLE_APPLICATION_CREDENTIALS=/config/gcloud/application_default_credentials.json \
    -e SNAKEMAKE_GOOGLEBATCH_PROJECT="${SNAKEMAKE_GOOGLEBATCH_PROJECT}" \
    -e GOOGLE_CLOUD_PROJECT="${GOOGLE_CLOUD_PROJECT}" \
    -e SNAKEMAKE_GOOGLEBATCH_REGION="${SNAKEMAKE_GOOGLEBATCH_REGION}" \
    snakemake-conda-apptainer:0.15-local \
    python3 /app/run.py -p local -c /app/genomes_download/config/config.yaml
