.. _`RNASeq`:

RNASEQ
======

Perform various analyses on RNASeq-bulk

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your directory containing VCF files

  cd /path/to/project

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/rnaseq/run.sh


Input/Output
------------


**Input:**

 
  
* Fastq files
  
 


**Output:**

 
  
* DGE + Fusions + reports
  
 









Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    ####################################
    ### General pipeline operations  ###
    ####################################


    include: "rules/000.common.smk"


    ##################
    ### Flag rules ###
    ##################


    onsuccess:
        shell("touch DONE && rm --force --verbose ON_GOING ERROR")
        logging.info("RNA-Seq bulk pipeline ended without any error.")


    onerror:
        shell("touch ERROR && rm --force --verbose ON_GOING DONE")
        logging.error("RNA-Seq bulk pipeline failed.")


    onstart:
        shell("touch ON_GOING && rm --force --verbose ERROR DONE")
        logging.info("Starting RNA-Seq bulk pipeline.")


    #############################
    ### Main target rule      ###
    ### Add new targets below ###
    #############################


    rule target:
        input:
            general_qc="data_output/multiqc/MultiQC.QC.html",
            quantification_qc="data_output/multiqc/MultiQC.Salmon.html",
            dge_qc=expand(
                "data_output/DEseq2/{comparison}/MultiQC.DEseq2.html",
                comparison=output_prefixes,
            ),
            chimera_qc="data_output/multiqc/MultiQC.Star.Chimera.html",


    ##########################
    ### Raw fastq handling ###
    ##########################

    """
    If you need to perform simple quality controls over
    your dataset, then use:
    --until quality_control_results
    """


    include: "rules/001.bigr_copy.smk"
    include: "rules/002.trimming.smk"
    include: "rules/003.fastq_screen.smk"


    rule quality_control_results:
        input:
            fastp=expand(
                "fastp/{ext}/{sample}.fastp.{ext}",
                sample=sample_list,
                ext=["html", "json"],
            ),
            fastq_screen=expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
                sample=sample_list,
                stream=streams,
                ext=["txt", "png"],
            ),
        output:
            report(
                "data_output/multiqc/MultiQC.QC.html",
                caption=str(workflow_source_dir / "reports" / "multiqc.rst"),
                category="Quality Controls",
            ),
        threads: 1
        resources:
            mem_mb=get_1gb_per_attempt,
            time_min=get_45min_per_attempt,
            tmpdir="tmp",
        retries: 2
        log:
            "logs/multiqc.salmon.log",
        wrapper:
            "bio/multiqc"


    #############################
    ### Salmon quantification ###
    #############################

    """
    If you need to run this pipeline and stop after a
    simple salmon quantification, use:
    --until salmon_quant_results
    """


    include: "rules/004.salmon.smk"


    """
    This snakefile calls salmon and aggregates counts in human-readable tables
    """


    # Call this fantom rule to access salmon quand only"
    rule salmon_quant_results:
        input:
            salmon_quant=expand(
                "salmon/pseudo_mapping/{sample}/quant.genes.sf", sample=sample_list
            ),
            tpm_table=expand(
                "data_output/Quantification/TPM.{feature}.tsv",
                feature=["transcripts", "genes"],
            ),
            raw_counts="data_output/Quantification/Raw.genes.tsv",
            fastp=expand(
                "fastp/{ext}/{sample}.fastp.{ext}",
                sample=sample_list,
                ext=["html", "json"],
            ),
            fastq_screen=expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
                sample=sample_list,
                stream=streams,
                ext=["txt", "png"],
            ),
        output:
            report(
                "data_output/multiqc/MultiQC.Salmon.html",
                caption=str(workflow_source_dir / "reports" / "multiqc.rst"),
                category="Quality Controls",
            ),
        threads: 1
        resources:
            mem_mb=get_1gb_per_attempt,
            time_min=get_45min_per_attempt,
            tmpdir="tmp",
        retries: 2
        log:
            "logs/multiqc.salmon.log",
        wrapper:
            "bio/multiqc"


    ####################
    ### Immunedeconv ###
    ####################

    """
    If you need to run this pipeline and stop after a
    simple immunedeconv, use:
    --until immunedeconv_results
    """


    include: "rules/005.immunedeconv.smk"


    rule immunedeconv_results:
        input:
            salmon_quant=lambda wildcards: expand(
                "salmon/pseudo_mapping/{sample}/quant.genes.sf",
                sample=sample_list,
            ),
            fastp=lambda wildcards: expand(
                "fastp/{ext}/{sample}.fastp.{ext}",
                sample=sample_list,
                ext=["html", "json"],
            ),
            fastq_screen=lambda wildcards: expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
                sample=sample_list,
                stream=streams,
                ext=["txt", "png"],
            ),
            config="immunedeconv/multiqc_config.yaml",
            graphs=expand(
                "data_output/{tool}/celltypes.{graph}.png",
                tool=config["immunedeconv"].get(
                    "tool_list",
                    [
                        "xcell",
                        "quantiseq",
                        "epic",
                        "mcpcounter",
                        "cibersort",
                        "cibersort_abs",
                    ],
                ),
                graph=["histogram", "dotplot"],
            ),
            tpm_table=expand(
                "data_output/Quantification/TPM.{feature}.tsv",
                feature=["transcripts", "genes"],
            ),
            raw_counts="data_output/Quantification/Raw.genes.tsv",
        output:
            protected("data_output/MultiQC/ImmuneDeconv.html"),
        threads: 1
        resources:
            mem_mb=get_2gb_per_attempt,
            time_min=get_45min_per_attempt,
            tmpdir="tmp",
        log:
            "logs/multiqc/immunedeconv.log",
        wrapper:
            "bio/multiqc"


    ####################
    ### Star chimera ###
    ####################

    """
    If you need to run this pipeline and stop after
    a simple star mapping for chimera and fusion reads,
    then use:
    --until star_fusion_results
    """


    include: "rules/011.mapping_qc.smk"
    include: "rules/012.star_chimera.smk"
    include: "rules/013.star_fusion.smk"
    include: "rules/019.rseqc.smk"


    rule star_fusion_results:
        input:
            mappings=expand(
                "star/{sample}/chimera/{sample}.{ext}",
                sample=sample_list,
                ext=["bam", "bam.bai"],
            ),
            #cram=expand("data_output/Mapping/chimera/{sample}.cram", sample=sample_list),
            cj=expand(
                "star/{sample}/chimera/{sample}.Chimeric.out.junction", sample=sample_list
            ),
            sj=expand("star/{sample}/chimera/{sample}.SJ.out.tab", sample=sample_list),
            log=expand("star/{sample}/chimera/{sample}.Log.out", sample=sample_list),
            logp=expand(
                "star/{sample}/chimera/{sample}.Log.progress.out", sample=sample_list
            ),
            logf=expand("star/{sample}/chimera/{sample}.Log.final.out", sample=sample_list),
            samtools_stats=expand(
                "samtools/stats/{sample}.chimera.stats", sample=sample_list
            ),
            fastp=expand(
                "fastp/{ext}/{sample}.fastp.{ext}",
                sample=sample_list,
                ext=["html", "json"],
            ),
            fastq_screen=expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
                sample=sample_list,
                stream=streams,
                ext=["txt", "png"],
            ),
            star_fusion=expand("star-fusions/{sample}/", sample=sample_list),
            read_distribution=expand(
                "rseqc/read_distribution/chimera/{sample}.read_distribution.tsv",
                sample=sample_list,
            ),
            tin=expand("rseqc/tin/chimera/{sample}.summary.txt", sample=sample_list),
            bam_stat=expand("rseqc/bam_stat/chimera/{sample}.txt", sample=sample_list),
            genebosycov=expand(
                "rseqc/gene_body_coverage/chimera/{sample}.geneBodyCoverage.{ext}",
                sample=sample_list,
                ext=["txt", "pdf"],
            ),
            txt=expand("rseqc/junction_annotation/chimera/{sample}.txt", sample=sample_list),
        threads: 1
        output:
            report(
                "data_output/multiqc/MultiQC.Star.Chimera.html",
                caption=str(workflow_source_dir / "reports" / "multiqc.rst"),
                category="Quality Controls",
            ),
        threads: 1
        resources:
            mem_mb=get_1gb_per_attempt,
            time_min=get_45min_per_attempt,
            tmpdir="tmp",
        retries: 2
        log:
            "logs/multiqc.star.chimera.log",
        wrapper:
            "bio/multiqc"


    ##############
    ### DESeq2 ###
    ##############

    """
    If you need to run this pipeline and stop after DESeq2,
    then use:
    --until deseq2_results
    """


    include: "rules/007.tximport.smk"
    include: "rules/008.deseq2.smk"
    include: "rules/009.deseq2_post_process.smk"


    rule deseq2_results:
        input:
            salmon_quant=lambda wildcards: expand(
                "salmon/pseudo_mapping/{sample}/quant.genes.sf",
                sample=samples_per_prefixes[wildcards.comparison],
            ),
            tpm_table=expand(
                "data_output/Quantification/TPM.{feature}.tsv",
                feature=["transcripts", "genes"],
            ),
            raw_counts="data_output/Quantification/Raw.genes.tsv",
            fastp=lambda wildcards: expand(
                "fastp/{ext}/{sample}.fastp.{ext}",
                sample=samples_per_prefixes[wildcards.comparison],
                ext=["html", "json"],
            ),
            fastq_screen=lambda wildcards: expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
                sample=samples_per_prefixes[wildcards.comparison],
                stream=streams,
                ext=["txt", "png"],
            ),
            deseq2_results=expand(
                "data_output/DEseq2/{comparison}/{content}_{comparison}.tsv",
                content=content_list,
                allow_missing=True,
            ),
            config="multiqc/{comparison}/multiqc_config.yaml",
            gene_plot=expand(
                "data_output/DEseq2/{comparison}/gene_plots/{gene}.png",
                gene=config.get("genes_of_interest", ["ENSG00000141510"]),
                allow_missing=True,
            ),
            plots=[
                "multiqc/{comparison}/pca_plot_mqc.png",
                "multiqc/{comparison}/volcanoplot_mqc.png",
                "multiqc/{comparison}/distro_expr_mqc.png",
                "multiqc/{comparison}/distro_mu_mqc.png",
                "multiqc/{comparison}/ma_plot_mqc.png",
                "multiqc/{comparison}/independent_filter_mqc.png",
                "multiqc/{comparison}/inde_theta_filter_mqc.png",
                "multiqc/{comparison}/pvalue_qc_mqc.png",
            ],
            deseq_rbt="data_output/DESeq2/{comparison}/Complete_html_table.tar.bz2",
        output:
            report(
                "data_output/DEseq2/{comparison}/MultiQC.DEseq2.html",
                caption=str(workflow_source_dir / "reports" / "multiqc.rst"),
                category="Quality Controls",
            ),
        threads: 1
        resources:
            mem_mb=get_1gb_per_attempt,
            time_min=get_45min_per_attempt,
            tmpdir="tmp",
        retries: 2
        log:
            "logs/multiqc.deseq2/{comparison}.log",
        wrapper:
            "bio/multiqc"


    #######################
    ### Variant Calling ###
    #######################


    include: "rules/010.star.smk"




Authors
-------


* Thibault Dayris
