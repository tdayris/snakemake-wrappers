# Somatic WES

1. [TLDR](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#tldr-run-this-pipeline)
1. [`design.tsv` file](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#design-file)
1. [`config.yaml` file](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#config-file)
1. [Classical use](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#classical-use)
1. [Quality controls](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#quality-controls)
1. [Mapping only](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#mapping-only)
1. [Variant calling](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#variant-calling)
1. [CNV calling](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#cnv-calling)
1. [Tumor Molecular Burden](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#tumor-molecular-burden)
1. [Micro satellites instability](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#miscosatellites-instability)
1. [Fusions](https://github.com/tdayris/snakemake-wrappers/tree/Unofficial/bigr_pipelines/wes_somatic#fusions) Under construction


# TLDR: run this pipeline

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input

# Run this pipeline
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh
```

# Design file

The design file contains paths to your samples. It should look like:

| Sample_id | Upstream_file_normal   | Downstream_file_normal   | Upstream_file_tumor      | Downstream_file_tumor   |
|-----------|------------------------|--------------------------|--------------------------|-------------------------|
| Sample1   | /path/file.1.fq.gz     | /path/file.2.fq.gz       | /path/file.1.fq.gz       | /path/file.2.fq.gz      |
| Sample2   | /path/file.1.fq.gz     | /path/file.2.fq.gz       | /path/file.1.fq.gz       | /path/file.2.fq.gz      |
| Sample3   | /path/file.1.fq.gz     | /path/file.2.fq.gz       | /path/file.1.fq.gz       | /path/file.2.fq.gz      |
| Sample4   | /path/file.1.fq.gz     | /path/file.2.fq.gz       | /path/file.1.fq.gz       | /path/file.2.fq.gz      |

Header **must** be the first line:

1. `Sample_id` **must** be in the header (first line). Column order does not matter.
1. `Upstream_file_normal` **must** be in the header (first line). Column order does not matter.
1. `Downstream_file_normal` **must** be in the header (first line). No single-ended analyzes here. Column order does not matter.
1. `Upstream_file_tumor` **must** be in the header (first line). Column order does not matter.
1. `Downstream_file_tumor` **must** be in the header (first line). No single-ended analyzes here. Column order does not matter.

No trio analysis available with this pipeline.

No germline analysis available with this pipeline.

No differential mutation analysis available with this pipeline.

# Config file

This file describes all optional parameters set in the pipeline.

Possible edition:

1. You may wan to select info fields available in final result table.

The complete list of info fields contains more than two hundred of possible annotation fields. Such a large table does not suit the weak computers belonging to our fellow biologists.

Any other parameter change may break the pipeline, do it at your own risks.

# Classical use

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input
ln -sfrv data_input/design.tsv . || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

# Run basic quality controls, keep intermediary files
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh QC --nt

# Choose wether to continue the analysis or not.

# Edit config file
# 1. with your annotations of interest
# 1. raise minimum mapping quality, etc
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh variants --nt

# Have a look at mapping and calling results

# Run CNV calling
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh cnv --nt

# Run TMB estimation
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh tnb --nt

# Run fusion analysis
# Under construction
```

# Quality controls 

## Pipeline

1. iRODS copy 
1. Fastp
1. FastqScreen
1. MultiQC

Use iRODS command to get your fastq files, then clean them with Fastp. Assess organism quality with FastqScreen, then aggregate quality reports with MultiQC.

## Command line

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input
ln -sfrv data_input/design.tsv design.tsv || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

# Run basic quality controls, keep intermediary files
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh QC --nt
```

# Mapping only

## Pipeline

1. iRODS copy (access iRODS collections and merge samples that were sequenced through multiple runs)
1. Fastp (trimm fastq reads)
1. bwa-mem2 (map against genome reference)
1. sambamba/samtools (fix mates, clean poor mapping qualities, remove duplicated, remove unmapped)
1. FastqScreen
1. MultiQC

## Command line

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input
ln -sfrv data_input/design.tsv design.tsv || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

# Run basic quality controls, keep intermediary files, map files
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh mapping --nt
```

# Variant Calling

## Pipeline

1. iRODS copy (access iRODS collections and merge samples that were sequenced through multiple runs)
1. Fastp (trim fastq reads)
1. bwa-mem2 (map against genome reference)
1. sambamba/samtools (fix mates, clean poor mapping qualities, remove duplicated, remove unmapped)
1. mutect2/GATK (call variants, assess multiple bias and go though HardFilters from GATK)
1. SnpEff/SnpSift/BCFtools/VCFtools (Annotated with 14 separate databases, mostly cancer oriented)
1. In-house scripts (Annotated with 4 databases built in Gustave Roussy or clinical data)
1. Facets (call CNV)
1. In-house scripts (Assess tumor molecular burden)
1. MSIsensor-pro (call MSI status)
1. FastqScreen
1. MultiQC

## Command line

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input
ln -sfrv data_input/design.tsv design.tsv || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

# Run basic quality controls, keep intermediary files, call variants, CNV, TMB, and MSI
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh variants --nt
```

# CNV calling

## Pipeline

1. iRODS copy (access iRODS collections and merge samples that were sequenced through multiple runs)
1. Fastp (trimm fastq reads)
1. bwa-mem2 (map against genome reference)
1. sambamba/samtools (fix mates, clean poor mapping qualities, remove duplicated, remove unmapped)
1. mutect2/GATK (call variants, assess multiple bias and go though HardFilters from GATK)
1. Facets (call CNV)
1. FastqScreen
1. MultiQC

## Command line

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input
ln -sfrv data_input/design.tsv design.tsv || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

# Run basic quality controls, keep intermediary files, call CNV
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh cnv --nt
```

# Tumor molecular burden

## Pipeline

1. iRODS copy (access iRODS collections and merge samples that were sequenced through multiple runs)
1. Fastp (trim fastq reads)
1. bwa-mem2 (map against genome reference)
1. sambamba/samtools (fix mates, clean poor mapping qualities, remove duplicated, remove unmapped)
1. mutect2/GATK (call variants, assess multiple bias and go though HardFilters from GATK)
1. In-house scripts (Assess tumor molecular burden)
1. FastqScreen
1. MultiQC

## Command line

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input
ln -sfrv data_input/design.tsv design.tsv || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

# Run basic quality controls, keep intermediary files, estimates TMB
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh tmb --nt
```

# Micro satellites instability

## Pipeline

1. iRODS copy (access iRODS collections and merge samples that were sequenced through multiple runs)
1. Fastp (trim fastq reads)
1. bwa-mem2 (map against genome reference)
1. sambamba/samtools (fix mates, clean poor mapping qualities, remove duplicated, remove unmapped)
1. mutect2/GATK (call variants, assess multiple bias and go though HardFilters from GATK)
1. MSIsensor-pro (call MSI status)
1. FastqScreen
1. MultiQC

## Command line

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input
ln -sfrv data_input/design.tsv design.tsv || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

# Run basic quality controls, keep intermediary files, estimate msi status
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh msi --nt
```

# Fusions

## Pipeline

1. iRODS copy (access iRODS collections and merge samples that were sequenced through multiple runs)
1. Fastp (trimm fastq reads)
1. STAR (mapping + fusions)
1. FusionAnnotator
1. FusionInspector
1. MultiQC

## Command line

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input
ln -sfrv data_input/design.tsv design.tsv || echo "No design found. Create it, or let the pipeline guess sample pairs (risky)"

# Run basic quality controls, keep intermediary files, search, check and annotate fusions
bash /mnt/beegfs/pipelines/unofficial-snakemake-wrappers/1.1.0/bigr_pipelines/wes_somatic/run.sh fusions --nt
```