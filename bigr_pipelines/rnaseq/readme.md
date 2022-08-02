# RNA-Seq bulk

This pipeline aims to perform classical RNASeq-bulk analyses:

1. TLDR
1. design.tsv file
1. config.yaml file
1. Classical use
1. Quality controls
1. Quantification
1. Differential Gene Expression
1. Aggregate factors
1. Ignore factors and/or levels
1. Perform only a subset of DGE
1. Fusion


Each subsection can be selected one by one, or the whole pipeline may be executed.

## TLDR: run this pipeline

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input

# Put design file in ${PWD}

# Run this pipeline
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh
```

## Design file

The design file contains a description of your samples and experimental design. It should look like:

| Sample_id | Upstream_file      | Downstream_file   | {ConditionA}   | {ConditionB} | ... |
+===========+====================+===================+================+==============+=====+
| Sample1   | /path/file.1.fq.gz |/path/file.2.fq.gz | Treated        | BatchA       | ... |
+-----------+--------------------+-------------------+----------------+--------------+-----+
| Sample1   | /path/file.1.fq.gz |/path/file.2.fq.gz | Treated        | BatchB       | ... |
+-----------+--------------------+-------------------+----------------+--------------+-----+
| Sample1   | /path/file.1.fq.gz |/path/file.2.fq.gz | Treated        | BatchA       | ... |
+-----------+--------------------+-------------------+----------------+--------------+-----+
| Sample1   | /path/file.1.fq.gz |/path/file.2.fq.gz | Untreated      | BatchB       | ... |
+-----------+--------------------+-------------------+----------------+--------------+-----+
| Sample1   | /path/file.1.fq.gz |/path/file.2.fq.gz | Untreated      | BatchA       | ... |
+-----------+--------------------+-------------------+----------------+--------------+-----+
| Sample1   | /path/file.1.fq.gz |/path/file.2.fq.gz | Untreated      | BatchB       | ... |
+-----------+--------------------+-------------------+----------------+--------------+-----+

`ConditionA` and `ConditionB` are exemple names. Put the name of your own conditions instead of ConditionA and ConditionB. Add as many conditions as you with below. Keep in mind that _all_ possible comparisons will be done (using `~{condition_name}`, see below for accounting batch effects).

Header line **must** be the first line:

1. `Sample_id` **must** be in the header (first line). Column order does not matter.
1. `Upstream_file` **must** be in the header (first line). Column order does not matter.
1. `Downstream_file` **must** be in the header (first line). No single-ended analyses here. Column order does not matter.
1. If `BatchEffect` is present in header (first line), and only if `BatchEffect` is present in header (first line), then the Differential Expression Analysis will take this column as batch effect (using `~BatchEffect+{condition_name}`). 

Only one batch effect can be considered. Do not call any of your conditions `BatchEffect` if you don't know what you are doing.

Each column will be considered as a factor. Each value in each cell will be considered as a level of its corresponding factor: 

1. A level with less than two replicates will be ignored. 
1. A factor with less than two levels will be ignored.


## Config.yaml

This file describes all optional parameters set in the pipeline. 

Possible edition:
1. You may want to modify the list of genes of interest by adding a list after the provided example. 
1. You may change the path to your design file by changing the default `design: design.tsv`.
1. You may want to analyse two factors as together (see: aggregate factors, TODO)
1. You may want to ignore factors present in your design file (see ignore factors, TODO)
1. You may want to run a subset of possible DGE (see subset DGE, TODO)

Do not change command line parameters if you do not know what you are doing. Do not set neither threading, temporary files, memory managment, nor IO files in command line. This has already been done by the pipeline.

Changing anything else shall be done at your own risks.


## Classical use

```{sh}
# Go to your project directory
cd /path/to/my/project/tmp

# Setup IO repositories
ln -sfrv ../data_output data_output || mkdir -pv data_output
ln -sfrv ../data_input data_input || mkdir -pv data_input

# Put design file in ${PWD}

# Run basic quality controls, keep intermediar files
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh QC --nt

# Choose wether to continue the analysis or not.

# Edit config file
# 1. with your genes of interest
# 1. subset list of DGE
# 1. aggregate factors, etc.
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh quant --nt

# Have a look at mapping rates and basic quantification QC

# Run Differential Gene Expression
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh dge --nt

# Have a look at the DGE results

# Run fusions
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh fusions --nt

# Have a look at the results

# Clean TEMPORARY files and temporary files only
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh --delete-temp-output
```

## Quality controls

### Pipeline

1. iRODS copy 
1. Fastp
1. FastqScreen
1. MultiQC

Use iRODS command to get your fastq files, then clean them with Fastp. Assess organism quality with FastqScreen, then aggregate quality reports with MultiQC.


### Command line

Command line argument order does not matter.

```{sh}
# GRCh38 / HG38
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh QC --nt
# GRCm38 / MM10
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh mm10 QC --nt
```

### Results

The only output is a MultiQC report. By default, all other output are deleted, use `--nt` to keep cleaned fastq files.

```
data_output/
├── multiqc
    ├── MultiQC.QC_data
    ├── multiqc_data.json
    │   ├── multiqc_fastp.txt
    │   ├── multiqc_fastq_screen.txt
    │   ├── multiqc_general_stats.txt
    │   ├── multiqc.log
    │   └── multiqc_sources.txt
    └── MultiQC.QC.html
```

## Quantification


### Pipeline

### Pipeline

1. iRODS copy 
1. Fastp
1. FastqScreen
1. Salmon
1. MultiQC


Use iRODS command to get your fastq files, then clean them with Fastp. Assess organism quality with FastqScreen. Salmon estimates transcripts abundance, then aggregate quality reports with MultiQC. Salmon results are annotated with a basic GTF.


### Command line

Command line argument order does not matter.

```{sh}
# GRCh38 / HG38
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh quant --nt
# GRCm38 / MM10
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh mm10 quant --nt
```

### Output

The repository `multiqc` contains two multiqc reports: basic quality controls and salmon. MultiQC.Salmon.html contains all information present in MultiQC.QC.html, and adds the results of Salmon.

The repository `Quantification` contains three quantification files:

1. The file: `Raw.genes.tsv` contains the raw gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is **not** normalized. Usefull for DESeq2, etc.
1. The file: `TPM.genes.tsv` contains TPM normalized gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is normalized by TPM. It does not take the factors and levels present in the `design.tsv` into account while normalizing.
1. The file: `TPM.transcripts.tsv` contains TPM normalized transcripts expression estimates, obtained Salmon. This is normalized by TPM. It does not take the factors and levels present in the `design.tsv` into account while normalizing.

```
data_output/
├── multiqc
│   ├── MultiQC.QC_data
│   │   ├── multiqc_data.json
│   │   ├── multiqc_fastp.txt
│   │   ├── multiqc_fastq_screen.txt
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc.log
│   │   └── multiqc_sources.txt
│   ├── MultiQC.QC.html
│   ├── MultiQC.Salmon_data
│   │   ├── multiqc_data.json
│   │   ├── multiqc_fastp.txt
│   │   ├── multiqc_fastq_screen.txt
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc.log
│   │   └── multiqc_sources.txt
│   └── MultiQC.Salmon.html
└── Quantification
    ├── Raw.genes.tsv
    ├── TPM.genes.tsv
    └── TPM.transcripts.tsv
```




## Differential Gene Expression

### Pipeline

1. iRODS copy 
1. Fastp
1. FastqScreen
1. Salmon
1. tximport
1. DESeq2
1. In-house scripts
1. MultiQC

### Command line

```{sh}
# GRCh38 / HG38
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh dge --nt
# GRCm38 / MM10
bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh mm10 dge --nt
```

### Results

The repository `multiqc` contains two multiqc reports: basic quality controls and salmon. MultiQC.Salmon.html contains all information present in MultiQC.QC.html, and adds the results of Salmon.

The repository `Quantification` contains three quantification files:

1. The file: `Raw.genes.tsv` contains the raw gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is **not** normalized. Usefull for DESeq2, etc.
1. The file: `TPM.genes.tsv` contains TPM normalized gene expression estimates, obtained from dummy summary of the expression of each transcripts counting for a same gene. This is normalized by TPM. It does not take the factors and levels present in the `design.tsv` into account while normalizing.
1. The file: `TPM.transcripts.tsv` contains TPM normalized transcripts expression estimates, obtained Salmon. This is normalized by TPM. It does not take the factors and levels present in the `design.tsv` into account while normalizing.


The `DEseq2` repository contains one sub-folder for each comparison between two levels, for each factor. Aggregation, sub-setting and factor ignore may reduce of increase the number of sub-forlders. See config file modification to know how to master the list of DESeq2 results. TODO

In each DESeq2 sub-folder:

1. The file `Complete_XXX.tsv` contains annotated results of DESeq2 with both differentially expressed genes under provided alpha threshold and the other ones. XXX being the name of the comparison.
1. The file `SortedOnLogFC_XXX.tsv` contains the differentially expressed genes only. They have been annotated and sorted on Log(FC). XXX being the name of the comparison.
1. The file `SortedOnPadj_XXX.tsv` contains the differentially expressed genes only. They have been annotated and sorted on adjusted P-Values. XXX being the name of the comparison.
1. The foler `gene_plots` contains a list of plots, one for each gene of interest. This is a highlight on its status, expression and per-condition weight.
1. The MultiQC report contains information about the samples belonging to the given comparison, and no other sample. It also contains Volcanoplot, PCAs, Expression weights, etc. Each of these plots is unique to this comparison, since expressions were normalized with DESeq2, which takes factors/levels into account.

```
data_output/
├── DEseq2
│   ├── DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2
│   │   ├── Complete_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
│   │   ├── gene_plots
│   │   │   └── ENSG00000141510.png
│   │   ├── MultiQC.DEseq2_data
│   │   │   ├── multiqc_data.json
│   │   │   ├── multiqc_fastp.txt
│   │   │   ├── multiqc_fastq_screen.txt
│   │   │   ├── multiqc_general_stats.txt
│   │   │   ├── multiqc.log
│   │   │   └── multiqc_sources.txt
│   │   ├── MultiQC.DEseq2.html
│   │   ├── SortedOnLogFC_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
│   │   └── SortedOnPadj_DGE_considering_factor_Batch_comparing_test_B1_vs_reference_B2.tsv
├── multiqc
│   ├── MultiQC.QC_data
│   │   ├── multiqc_data.json
│   │   ├── multiqc_fastp.txt
│   │   ├── multiqc_fastq_screen.txt
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc.log
│   │   └── multiqc_sources.txt
│   ├── MultiQC.QC.html
│   ├── MultiQC.Salmon_data
│   │   ├── multiqc_data.json
│   │   ├── multiqc_fastp.txt
│   │   ├── multiqc_fastq_screen.txt
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc.log
│   │   └── multiqc_sources.txt
│   └── MultiQC.Salmon.html
└── Quantification
    ├── Raw.genes.tsv
    ├── TPM.genes.tsv
    └── TPM.transcripts.tsv
```