.. _`gatk_bqsr`:

GATK_BQSR
=========

Compute and apply base recalibration on BAM files with GATK


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_gatk_bqsr = {
        "threads": 20,
        "genome": "reference/genome.fa",
        "dbsnp": "reference/dbsnp.vcf.gz"
    }

    def get_fasta_index_from_genome_path(genome_path: str) -> str:
        return genome_path + ".fai"

    def get_fasta_dict_from_genome_path(genome_path: str) -> str:
        return ".".join(genome_path.split(".")[:-1]) + ".dict"


    def get_vcf_tbi_from_db_path(db_path: str) -> str:
        return db_path + ".tbi"

    try:
        if config == dict():
            config = default_config_gatk_bqsr
    except NameError:
        config = default_config_gatk_bqsr

    """
    This rule applies the BQSR to the mapped reads
    """
    rule gatk_apply_baserecalibrator:
        input:
            bam="mapped/{sample}.bam",
            bam_index="mapped/{sample}.bam.bai",
            ref=config["genome"],
            ref_idx=get_fasta_index_from_genome_path(config["genome"]),
            ref_dict=get_fasta_dict_from_genome_path(config["genome"]),
            recal_table="gatk/recal_data_table/{sample}.grp"
        output:
            bam="gatk/recal_bam/{sample}.bam"
        message: "Applying BQSR on {wildcards.sample} with GATK"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4096, 20240),
            time_min=lambda wildcards, attempt: attempt * 60
        log:
            "logs/gatk/applybqsr/{sample}.log"
        params:
            extra=""
        wrapper:
            "0.72.0-560-g28998a654/bio/gatk/applybqsr"


    """
    This rule computes BQSR on mapped reads, given a knoledge database
    """
    rule gatk_compute_baserecalibration_table:
        input:
            bam="mapped/{sample}.bam",
            bam_index="mapped/{sample}.bam.bai",
            ref=config["genome"],
            ref_idx=get_fasta_index_from_genome_path(config["genome"]),
            ref_dict=get_fasta_dict_from_genome_path(config["genome"]),
            known=config["dbsnp"],
            known_idx=get_vcf_tbi_from_db_path(config["dbsnp"])
        output:
            recal_table=temp("gatk/recal_data_table/{sample}.grp")
        message: "Compute BQSR table from {wildcards.sample} with GATK"
        threads: config.get("threads", 20)
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4048, 15360),
            time_min=lambda wildcards, attempt: attempt * 120
        log:
            "logs/gatk3/compute_bqsr/{sample}.log"
        params:
            extra=""
        wrapper:
            "0.72.0-560-g28998a654/bio/gatk/baserecalibrator"


    """
    This rule indexes the input genome sequence with Samtools. It is not
    explicitely requested by GATK, but it will crash if the genome sequence
    is not indexed.

    This rule is cached since it should be used only once per reference sequence
    """
    rule samtools_faidx:
        input:
            config["genome"]
        output:
            get_fasta_dict_from_genome_path(config["genome"])
        message: "Indexing reference fasta with Samtools"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 4098),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            ""
        log:
            "logs/samtools/faidx/{genome}.log"
        wrapper:
            "0.72.0-560-g28998a654/bio/samtools/faidx"


    """
    This rule creates a sequence dictionnary from a genome sequnece. It is not
    explicitely requested by GATK, but it will crash if the genome sequence
    is not indexed.

    This rule is cached since it should be used only once per reference sequence
    """
    rule picard_create_sequence_dictionnary:
        input:
            config["genome"]
        output:
            get_fasta_dict_from_genome_path(config["genome"])
        message: "Creating sequence dictionnary over reference genome with Picard"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 2048, 8192),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            ""
        log:
            "logs/picard/create_sequence_dictionnary/{genome}.log"
        wrapper:
            "0.72.0-560-g28998a654/bio/picard/createsequencedictionary"


    """
    This rule creates a TBI index for the known VCF file. It is not
    explicitely requested by GATK, but it will crash if the genome sequence
    is not indexed.

    This rule is cached since it should be used only once per reference sequence
    """
    rule tabix_index:
        input:
            config["dbsnp"]
        output:
            get_vcf_tbi_from_db_path(config["dbsnp"])
        message: "Indexing kown variants with Tabix"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 10240),
            time_min=lambda wildcards, attempt: attempt * 60
        params:
            "-p vcf"
        log:
            "logs/tabix/index/{known}.log"
        wrapper:
            "0.72.0-560-g28998a654/bio/tabix"


    """
    This rule indexes the bam file with Samtools. It is not
    explicitely requested by GATK, but it will crash if the genome sequence
    is not indexed.
    """
    rule samtools_index:
        input:
            "mapped/{sample}.bam"
        output:
            "mapped/{sample}.bam.bai"
        message: "Indexing mapped reads of {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=1536,
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/samtools/sort/{sample}.log"
        wrapper:
            "0.72.0-560-g28998a654/bio/samtools/index"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/gatk/baserecalibrator`

* :ref:`bio/gatk/applybqsr`

* :ref:`bio/picard/createsequencedictionary`

* :ref:`bio/samtools/faidx`

* :ref:`bio/samtools/index`

* :ref:`bio/tabix`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

From: https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-

> You should almost always perform recalibration on your sequencing data. In human data, given the exhaustive databases of variation we have available, almost all of the remaining mismatches -- even in cancer -- will be errors, so it's super easy to ascertain an accurate error model for your data, which is essential for downstream analysis. For non-human data it can be a little bit more work since you may need to bootstrap your own set of variants if there are no such resources already available for you organism, but it's worth it.

Warning:

* Bam files must have read groups




Authors
-------


* Thibault Dayris

