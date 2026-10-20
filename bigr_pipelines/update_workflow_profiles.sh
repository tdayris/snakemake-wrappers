#!/usr/bin/bash
set -euo 'pipefail'
shopt -s 'nullglob'

echo "Updating genomes_download pipeline"
echo " L Updating gcloud workflow-profile"
mkdir -pv genomes_download/workflow-profiles/{gcloud,slurm}
cat '../meta/bio/genome_statistics/resources/workflow-profiles/gcloud/config.yaml' \
    > 'genomes_download/workflow-profiles/gcloud/config.yaml'
sed '1d' '../meta/bio/download_genome_sequences/resources/workflow-profile/gcloud/config.yaml' \
    >> 'genomes_download/workflow-profiles/gcloud/config.yaml'
sed '1d' '../meta/bio/download_genome_annotation/resources/workflow-profile/gcloud/config.yaml' \
    >> 'genomes_download/workflow-profiles/gcloud/config.yaml'

echo " L Updating slurm workflow-profile"
cat '../meta/bio/genome_statistics/resources/workflow-profiles/slurm/config.yaml' \
    > 'genomes_download/workflow-profiles/slurm/config.yaml'
sed '1d'    '../meta/bio/download_genome_sequences/resources/workflow-profile/slurm/config.yaml' \
    >> 'genomes_download/workflow-profiles/slurm/config.yaml'
sed '1d'    '../meta/bio/download_genome_annotation/resources/workflow-profile/slurm/config.yaml' \
    >> 'genomes_download/workflow-profiles/slurm/config.yaml'

echo "Process over"
