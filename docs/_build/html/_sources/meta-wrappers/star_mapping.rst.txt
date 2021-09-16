.. _`star_mapping`:

STAR_MAPPING
============

Map your reads with star


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_star_config = {
        "fasta": "/path/to/fasta",
        "gtf": "/path/to/gtf"
    }


    try:
        if config == dict():
            config = default_star_config
    except NameError:
        config = default_star_config


    ########################
    ### Quality controls ###
    ########################
    rule picard_metrics:
        input:
            ref=config["fasta"],
            ref_dict=get_dict(config["fasta"]),
            ref_fai=get_fai(config["fasta"]),
            bam="bam/star/{sample}.bam"
            bai=get_bai("bam/star/{sample}.bam")
        output:
            "picard/alignment_summary/{sample}.summary.txt"
        message:
            "Collecting alignment summary for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/picard/alignment_summary/{sample}.log"
        wrapper:
            "bio/picard/collectalignmentsummarymetrics"


    rule samtools_view:
        input:
            "star/{sample}/Aligned.out.sam"
        output:
            "bam/star/{sample}.bam"
        message:
            "Filtering, sorting and compressing {wildcards.sample}"
        threads: 10
        resources:
            mem_mb=lambda: wildcards, threads: threads * 720,
            time_min=lambda: wildcards, attempt: attempt * 30,
            tmpdir="tmp"
        params:
            extra="-bhq 10"
        log:
            "logs/samtools/view/{sample}.log"
        wrapper:
            "bio/samtools/view"


    ###################
    ### STAR itself ###
    ###################


    rule star_mapping:
        input:
            fq1="reads/{sample}.1.fq.gz",
            fq2="reads/{sample}.2.fq.gz",
            index="star/index"
        output:
            pipe("star/{sample}/Aligned.out.sam"),
            tmp = temp(directory("star/tmp/{sample}.tmp"))
        message:
            "Mapping {wildcards.sample} to reference genome with STAR"
        threads: 20
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 5210 + 40960,
            time_min=lambda wildcards, attempt: attempt * 60,
            tmpdir="tmp"
        params:
            extra = (
                " --outFilterType BySJout "
                " --outFilterMultimapNmax 20 "
                " --alignSJoverhangMin 8 "
                " --alignSJDBoverhangMin 1 "
                " --outFilterMismatchNmax 999 "
                " --outFilterMismatchNoverReadLmax 0.04 "
                " --alignIntronMin 20 "
                " --alignIntronMax 1000000 "
                " --alignMatesGapMax 1000000 "
                " --twopassMode Basic "
                " --outTmpDir star/tmp/{sample}.tmp/ "
                " --outSAMattributes All "
            )
        log:
            "logs/star/align/{sample}.log"
        wrapper:
            "bio/star/align"


    rule star_index:
        input:
            fasta = config["fasta"]
            gtf = config["gtf"]
        output:
            directory("star/index")
        message:
            "Indexing genome with STAR"
        threads: 20
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 5120 + 40960,
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir="tmp"
        params:
            extra="",
            sjdbOverhang="100"
        log:
            "logs/star/index.log"
        wrapper:
            "bio/star/index"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/star/index`

* :ref:`bio/star/align`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

