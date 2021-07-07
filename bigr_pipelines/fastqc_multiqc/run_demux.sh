#!/usr/bin/env bash

bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/fastqc_multiqc/run.sh -p demux output/multiqc.html "$@"
