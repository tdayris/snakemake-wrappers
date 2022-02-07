.. _`celltype_deconvolution_rnaseq (under development)`:

CELLTYPE_DECONVOLUTION_RNASEQ (UNDER DEVELOPMENT)
=================================================

Perform trimming and quantification on RNASeq, then cell type opulation deconvolution

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/salmon_quant/run.sh

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/salmon_quant/run.sh hg38


Input/Output
------------


**Input:**

 
  
* Fastq files
  
 
  
* Fasta-formatted Genome sequence
  
 
  
* Fasta-formatted transcriptome sequence
  
 
  
* GTF formatted genome annotation
  
 


**Output:**

 
  
* Both PNG and TSV formatted deconvolution results
  
 
  
* Salmon quantification
  
 
  
* Quality controls
  
 
  
* Trimmed fastq files
  
 







Notes
-----

Prerequisites:

* A TSV formatted design file, *named 'design.tsv'* with the following columns:

.. list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - Sample_id
    - Upstream_fastq
    - Downstream_fastq
  * - Name of the Sample1
    - Path to upstream fastq file
    - Path to downstream fastq file
  * - Name of the Sample2
    - Path to upstream fastq file
    - Path to downstream fastq file
  * - ...
    - ...
    - ...





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    from snakemake.utils import min_version
    from pathlib import Path
    from yaml import dump
    min_version("6.0")

    import sys

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../common/python")
    sys.path.append(common)

    from file_manager import *
    from files_linker import *
    from write_yaml import *
    from messages import message

    logging.basicConfig(
        filename="snakemake.celltype_deconvolution_rnaseq.log",
        filemode="w",
        level=logging.DEBUG
    )

    default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
    configfile: get_config(default_config)
    design = get_design(os.getcwd(), search_fastq_pairs)
    #print(design)

    if "Sample_id" in design.columns.tolist():
        design.set_index("Sample_id", inplace=True)

    fastq_links = link_fq(
        design.index,
        design.Upstream_file,
        design.Downstream_file
    )

    localrules: get_cibersort


    rule target:
        input:
            graphs = expand(
                "{tool}/celltypes.{plot}.png",
                tool=["cibersort_abs", "cibersort", "quantiseq", "xcell", "epic"],
                plot=["hist", "dotplot"]
            ),
            text = expand(
                "{tool}/celltypes.{ext}",
                tool=["cibersort_abs", "cibersort", "quantiseq", "xcell", "epic"],
                ext=["tsv", "RDS"]
            ),
            plotdirs = expand(
                "{tool}/celltypes.dotplots",
                tool=["cibersort_abs", "cibersort", "quantiseq", "xcell", "epic"]
            )

    ################################
    ### Deconvolution: CIBERSORT ###
    ################################


    rule cibersort_abs:
        input:
            expr_mat="immunedeconv/TPM.tsv",
            cibersort_binary="CIBERSORT.R",
            cibersort_mat="LM22.txt"
        output:
            histogram="cibersort_abs/celltypes.hist.png",
            dotplot="cibersort_abs/celltypes.dotplot.png",
            tsv="cibersort_abs/celltypes.tsv",
            rds="cibersort_abs/celltypes.RDS",
            plotdir=directory("cibersort_abs/celltypes.dotplots")
        message:
            "Using Cibersort-absolute to deconvolute expression into cell types"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmp="tmp"
        params:
            gene_col="Hugo_ID"
        log:
            "logs/immunedeconv/cibersort_abs.log"
        wrapper:
            "bio/immunedeconv/cibersort-abs"


    rule cibersort:
        input:
            expr_mat="immunedeconv/TPM.tsv",
            cibersort_binary="CIBERSORT.R",
            cibersort_mat="LM22.txt"
        output:
            histogram="cibersort/celltypes.hist.png",
            dotplot="cibersort/celltypes.dotplot.png",
            tsv="cibersort/celltypes.tsv",
            rds="cibersort/celltypes.RDS",
            plotdir=directory("cibersort/celltypes.dotplots")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmp="tmp"
        message:
            "Using cibersort to deconvolute expression into cell types"
        log:
            "logs/immunedeconv/Cibersort.log"
        params:
            gene_col="Hugo_ID"
        wrapper:
            "bio/immunedeconv/cibersort"


    rule get_cibersort:
        input:
            cibersort_binary="/mnt/beegfs/software/cibersort/1.0.6/CIBERSORT.R",
            cibersort_mat="/mnt/beegfs/software/cibersort/1.0.6/LM22.txt"
        output:
            cibersort_binary=temp("CIBERSORT.R"),
            cibersort_mat=temp("LM22.txt")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        message:
            "Gathering Cibersort requirements"
        log:
            "logs/cibersort.bin.log"
        params:
            rsync = "--verbose --checksum --human-readable",
            chmod = "u+x"
        shell:
            "rsync {params.rsync} {input.cibersort_binary} {output.cibersort_binary} > {log} 2>&1 && "
            "chmod {params.chmod} {output.cibersort_binary} >> {log} 2>&1 && "
            "rsync {params.rsync} {input.cibersort_mat} {output.cibersort_mat} >> {log} 2>&1"


    #############################
    ### Deconvolution: OTHERS ###
    #############################


    rule mcpcounter:
        input:
            expr_mat="immunedeconv/TPM.tsv"
        output:
            histogram="mcpcounter/celltypes.hist.png",
            dotplot="mcpcounter/celltypes.dotplot.png",
            tsv="mcpcounter/celltypes.tsv",
            rds="mcpcounter/celltypes.RDS",
            plotdir=directory("mcpcounter/celltypes.dotplots")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmp="tmp"
        message:
            "Using MCP-Counter to deconvolute expression into cell types"
        params:
            gene_col="Hugo_ID"
        log:
            "logs/immunedeconv/mcpcounter.log"
        wrapper:
            "bio/immunedeconv/mcpcounter"


    rule epic:
        input:
            expr_mat="immunedeconv/TPM.tsv"
        output:
            histogram="epic/celltypes.hist.png",
            dotplot="epic/celltypes.dotplot.png",
            tsv="epic/celltypes.tsv",
            rds="epic/celltypes.RDS",
            plotdir=directory("epic/celltypes.dotplots")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmp="tmp"
        params:
            gene_col="Hugo_ID"
        message:
            "Using EPIC to deconvolute expression into cell types"
        log:
            "logs/immunedeconv/epic.log"
        wrapper:
            "bio/immunedeconv/epic"


    rule quantiseq:
        input:
            expr_mat="immunedeconv/TPM.tsv"
        output:
            histogram="quantiseq/celltypes.hist.png",
            dotplot="quantiseq/celltypes.dotplot.png",
            tsv="quantiseq/celltypes.tsv",
            rds="quantiseq/celltypes.RDS",
            plotdir=directory("quantiseq/celltypes.dotplots")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmp="tmp"
        message:
            "Using QuantiSeq to deconvolute expression into cell types"
        params:
            gene_col="Hugo_ID"
        log:
            "logs/immunedeconv/quantiseq.log"
        wrapper:
            "bio/immunedeconv/quantiseq"


    rule xcell:
        input:
            expr_mat="immunedeconv/TPM.tsv"
        output:
            histogram="xcell/celltypes.hist.png",
            dotplot="xcell/celltypes.dotplot.png",
            tsv="xcell/celltypes.tsv",
            rds="xcell/celltypes.RDS",
            plotdir=directory("xcell/celltypes.dotplots")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmp="tmp"
        message:
            "Using xCell to deconvolute expression into cell types"
        params:
            gene_col="Hugo_ID"
        log:
            "logs/immunedeconv/xcell.log"
        wrapper:
            "bio/immunedeconv/xcell"

    ######################
    ### Quantification ###
    ######################


    rule subset_gene_counts:
        input:
            table="salmon_quant/salmon/TPM.genes.tsv"
        output:
            table="immunedeconv/TPM.tsv"
        message:
            "Formatting counts for ImmuneDeconv"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 512,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/immunedeconv/filter_gene_counts.log"
        params:
            drop_column = ["target_id"],
            set_index = "Hugo_ID",
            drop_duplicated_lines = True,
            keep_index = True,
            drop_duplicated_index = True,
            override_previous_index = True,
            filters = [
                ["Hugo_ID", "!=", "Unknown"]
            ]
        wrapper:
            "bio/pandas/filter_table"


    rule run_salmon_pipeline:
        input:
            config = "config.yaml"
        output:
            qc = "salmon_quant/multiqc/MultiQC.html",
            gene_counts = "salmon_quant/salmon/TPM.genes.tsv",
            transcript_counts = "salmon_quant/salmon/TPM.transcripts.tsv",
            results = directory("salmon_quant/results_to_upload")
        message:
            "Quantifying reads over the genome"
        handover: True
        priority: 10
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 60 + 45,
            tmpdir="tmp"
        log:
            "logs/salmon_quant_pipeline.log"
        params:
            ln = "--force --symbolic --relative",
            mk = "--verbose --parents"
        shell:
            "mkdir {params.mk} salmon_quant && "
            "ln {params.ln} design.tsv salmon_quant/design.tsv && "
            "cd salmon_quant/ && "
            "bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/salmon_quant/run.sh"




Authors
-------


* Thibault Dayris
