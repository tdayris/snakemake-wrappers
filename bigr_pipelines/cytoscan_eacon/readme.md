.. cytoscan_eacon:

# Cytoscan EaCoN

## Summary

This pipeline aims to perform classical CEL file analysis on Cytoscan array. It produces a HTML report, a set of images, and a Genomic Instability Score (GIS).

1. [TLDR]()
1. [design.tsv file]()
1. [config.yaml file]()
1. [Classical use]()
1. [ASCN/GIS rules on error]()

 ## TLDR: run this pipeline

 ```{sh}
 ## Enter batch name
 BATCH_NAME="MyBatch"

 ## Go to your project
 mkdir -vp "/mnt/beegfs/scratch/bioinfo_core/Oncoscans_Cytoscans/${BATCH_NAME}"
 cd "/mnt/beegfs/scratch/bioinfo_core/Oncoscans_Cytoscans/${BATCH_NAME}"

 ## Copy CEL files from colibri
 rsync -cvrhP ${USER}@colibri:/data/share/analyses/AffymetrixIGR/GenetiqueDesTumeurs/${BATCH_NAME} .

 ## Run this pipeline
 bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.1/bigr_pipelines/cytoscan_eacon/run.sh

 # Send data to NextCloud
 module load java/1.8.0_281-jdk

# WARNING: Create destination dir on Nextcloud
java -jar /mnt/beegfs/software/agent-nextcloud/agent-nextcloud-1.0.2.jar UPLOAD ${PWD} Onco_Cyto/${BATCH_NAME} $USER
```