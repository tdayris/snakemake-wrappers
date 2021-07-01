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

    fastq_links = link_fq(
        design.Sample_id,
        design.Upstream_file,
        design.Downstream_file
    )

    ruleorder: salmon_meta_salmon_quant_paired > salmon_quant_paired


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
                plot=["tsv", "rds"]
            ),
            plotdirs = expand(
                "{tool}/celltypes.dotplots",
                tool=["cibersort_abs", "cibersort", "quantiseq", "xcell", "epic"]
            )

    #####################
    ### Deconvolution ###
    #####################


    rule cibersort_abs:
        input:
            quant="immunedeconv/TPM.tsv"
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
            time_min=lambda wildcards, attempt: attempt * 15
        params:
            cibersort_binary="/path/to/cibersort/non/free",
            cibersort_mat="/path/to/celltypes.mat"
        log:
            "logs/immunedeconv/cibersort_abs.log"
        wrapper:
            "/bio/immunedeconv/cibersort-abs"


    use rule cibersort_abs as cibersort with:
        output:
            histogram="cibersort/celltypes.hist.png",
            dotplot="cibersort/celltypes.dotplot.png",
            tsv="cibersort/celltypes.tsv",
            rds="cibersort/celltypes.RDS",
            plotdir=directory("cibersort/celltypes.dotplots")
        message:
            "Using cibersort to deconvolute expression into cell types"
        log:
            "logs/immunedeconv/Cibersort.log"
        wrapper:
            "/bio/immunedeconv/cibersort"


    use rule cibersort_abs as mcpcounter with:
        output:
            histogram="mcpcounter/celltypes.hist.png",
            dotplot="mcpcounter/celltypes.dotplot.png",
            tsv="mcpcounter/celltypes.tsv",
            rds="mcpcounter/celltypes.RDS",
            plotdir=directory("mcpcounter/celltypes.dotplots")
        message:
            "Using MCP-Counter to deconvolute expression into cell types"
        log:
            "logs/immunedeconv/mcpcounter.log"
        wrapper:
            "/bio/immunedeconv/mcpcounter"


    use rule cibersort_abs as epic with:
        output:
            histogram="epic/celltypes.hist.png",
            dotplot="epic/celltypes.dotplot.png",
            tsv="epic/celltypes.tsv",
            rds="epic/celltypes.RDS",
            plotdir=directory("epic/celltypes.dotplots")
        message:
            "Using EPIC to deconvolute expression into cell types"
        log:
            "logs/immunedeconv/epic.log"
        wrapper:
            "/bio/immunedeconv/epic"


    use rule cibersort_abs as quantiseq with:
        output:
            histogram="quantiseq/celltypes.hist.png",
            dotplot="quantiseq/celltypes.dotplot.png",
            tsv="quantiseq/celltypes.tsv",
            rds="quantiseq/celltypes.RDS",
            plotdir=directory("quantiseq/celltypes.dotplots")
        message:
            "Using QuantiSeq to deconvolute expression into cell types"
        log:
            "logs/immunedeconv/quantiseq.log"
        wrapper:
            "/bio/immunedeconv/quantiseq"


    use rule cibersort_abs as xcell with:
        output:
            histogram="xcell/celltypes.hist.png",
            dotplot="xcell/celltypes.dotplot.png",
            tsv="xcell/celltypes.tsv",
            rds="xcell/celltypes.RDS",
            plotdir=directory("xcell/celltypes.dotplots")
        message:
            "Using xCell to deconvolute expression into cell types"
        log:
            "logs/immunedeconv/xcell.log"
        wrapper:
            "/bio/immunedeconv/xcell"


    #########################
    ### Extracting counts ###
    #########################

    rule gene_counts:
        input:
            table="salmon/TPM.genes.tsv"
        output:
            table="immunedeconv/TPM.tsv"
        message: "Cleaning salmon quantification file"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 5
        logs:
            "logs/pandas/cleaning.log"
        params:
            drop_column=["Strand", "Start", "Stop"],
            drop_null_lines=True,
            drop_na_lines=True,
            drop_duplicated_lines=True
        wrapper:
            "bio/pandas/filter_table"


    rule pandas_merge_salmon_tr:
        input:
            quant = expand(
                "salmon/pseudo_mapping/{sample}/quant.sf",
                sample=design.index.tolist()
            ),
            tx2gene = "tximport/transcripts2genes.tsv"
        output:
            tsv = "salmon/TPM.{content}.tsv"
        message:
            "Testing pandas merge salmon"
        params:
            header = True,
            position = True,
            gencode = False,
            genes = lambda wildcards: wildcards.content == "genes"
        log:
            "logs/pandas_merge_salmon/{content}.log"
        wrapper:
            "/bio/pandas/salmon"

    #############################
    ### Salmon quantification ###
    #############################

    salmon_config = {
        "genome": config["ref"]["genome"],
        "transcriptome": config["ref"]["transcriptome"],
        "gtf": config["ref"]["gtf"],
        "salmon_libtype": config["params"]["salmon_libtype"],
        "salmon_quant_extra": config["params"]["salmon_quant_extra"],
        "salmon_index_extra": config["params"]["salmon_index_extra"]
    }


    module salmon_meta:
        snakefile: "../../meta/bio/salmon/test/Snakefile"
        config: salmon_config


    rule multiqc:
        input:
            salmon=expand(
                "salmon/pseudo_mapping/{sample}/quant.sf",
                sample=design["Sample_id"]
            ),
            html=expand(
                "fastp/html/pe/{sample}.fastp.html",
                sample=design["Sample_id"]
            ),
            json=expand(
                "fastp/json/pe/{sample}.fastp.json",
                sample=design["Sample_id"]
            )
        output:
            report(
                "multiqc/MultiQC.html",
                caption="../common/reports/multiqc.rst",
                category="Quality Controls"
            )
        message:
            "Aggregating quality reports from Fastp and Salmon"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        log:
            "logs/multiqc.log"
        wrapper:
            "/bio/multiqc"


    use rule * from salmon_meta as salmon_meta_*


    use rule salmon_quant_paired from salmon_meta with:
        output:
            quant=report(
                "salmon/pseudo_mapping/{sample}/quant.sf",
                category="2. Raw Salmon output",
                caption="../../common/reports/salmon_quant.rst"
            ),
            lib="salmon/pseudo_mapping/{sample}/lib_format_counts.json",
            mapping=temp("salmon/bams/{sample}.bam")


    ############################
    ### FASTP FASTQ CLEANING ###
    ############################

    rule fastp_clean:
        input:
            sample=expand(
                "reads/{sample}.{stream}.fq.gz",
                stream=["1", "2"],
                allow_missing=True
            ),
        output:
            trimmed=expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
            html="fastp/html/pe/{sample}.fastp.html",
            json=temp("fastp/json/pe/{sample}.fastp.json")
        message: "Cleaning {wildcards.sample} with Fastp"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            adapters=config["params"].get("fastp_adapters", None),
            extra=config["params"].get("fastp_extra", "")
        log:
            "logs/fastp/{sample}.log"
        wrapper:
            "/bio/fastp"


    #################################################
    ### Gather files from iRODS or mounting point ###
    #################################################

    rule bigr_copy:
        output:
            "reads/{sample}.{stream}.fq.gz"
        message:
            "Gathering {wildcards.sample} fastq file ({wildcards.stream})"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            input=lambda wildcards, output: fastq_links[output[0]]
        log:
            "logs/bigr_copy/{sample}.{stream}.log"
        wrapper:
            "/bio/BiGR/copy"




Authors
-------


* Thibault Dayris
