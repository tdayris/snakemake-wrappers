These pipelines are supposed to work on BiGR Flamingo only.

# Installation and further updates:
WARNING: You do not have to install anything. You do not need to update anything. It has already been done for you.
*   `bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/common/install.sh"`

# Cytoscan:
1.  Go to your repository with CEL Cytoscan files: `cd /path/to/CEL/dir`
2.  Run pipeline: `bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/cytoscan_eacon/run.sh"`

->Test: `bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/cytoscan_eacon/test.sh"`

# Oncoscan:
1.  Go to your repository with CEL Oncoscan files: `cd /path/to/CEL/dir`
2.  Run pipeline: `bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/oncoscan_eacon/run.sh"`

->Test: `bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/oncoscan_eacon/test.sh"`

# Variant Calling Ampliseq (HG19):
1. Go to your repository with Fastq files: `cd /path/to/fastq/dir`
2. Run pipeline: `bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_ampliseq/run.sh"`

->Test: `bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_ampliseq/test.sh"`

# Variant Calling Ampliseq (HG38):
1. Go to your repository with Fastq files: `cd /path/to/fastq/dir`
2. Run pipeline: `bash "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_ampliseq/run.sh" hg38`

Warning: Ampliseq capture kits are made based on HG19.
