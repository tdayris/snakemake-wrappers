#!/usr/bin/env bash

set -e

#TODO: Read pirestParameters.tsv
# -> Get datasetIdBED col, ideally col 2
# -> Get genomeVersion col for genome version, ideally col #
# -> Set bed / genome versions
# -> Strip quotes
#
# Bientôt: GRCh39 #Pas publié
#
#Pour le moment:
#- GRCh38 # OK
#- GRCh19 # En cours de test
#- GRCm9  # Pas fait
#- GRCm38 # OK

# Col1 = Sample_id
# Col2 = datasetIdBED
# Col3 = Upstream_file_normal
# Col4 = Downstream_file_normal
# Col5 = Upstream_file_tumor
# Col6 = Downstream_file_tumor
# Col7 = genomeVersion

bed_file=$(cut -f2 "pirestParameters.tsv" | sort | uniq)
genome_version=$(cut -f7 "pirestParameters.tsv" | sort | uniq)

if [ "${#bed_file[@]}" -gt 1 ]; then
  echo "There was more than one single BED required for this analysis"
  #exit 1
fi

if [ "${#genome_version[@]}" -gt 1 ]; then
  echo "There was more than one genome version provided for this analysis"
  #exit 2
fi

if [ "${genome_version}" == "GRCh38" ]; then
  echo bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_somatic/run.sh hg38
elif [ "${genome_version}" == "GRCh37" ]; then
  echo bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_somatic/run.sh hg19
elif [ "${genome_version}" == "GRCm38" ]; then
  echo bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_somatic/run.sh mm10
else
  echo "Unknown genome version"
fi
