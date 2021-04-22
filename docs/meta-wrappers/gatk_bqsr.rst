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
        "genome": "reference/genome.fasta",
        "dbsnp": "reference/dbsnp.vcf.gz"
    }

    def get_fai(genome_path: str) -> str:
        return genome_path + ".fai"

    def get_dict(genome_path: str) -> str:
        return ".".join(genome_path.split(".")[:-1]) + ".dict"


    def get_tbi(db_path: str) -> str:
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
            ref_idx=get_fai(config["genome"]),
            ref_dict=get_dict(config["genome"]),
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
            "/bio/gatk/applybqsr"


    """
    This rule computes BQSR on mapped reads, given a knoledge database
    """
    rule gatk_compute_baserecalibration_table:
        input:
            bam="mapped/{sample}.bam",
            bam_index="mapped/{sample}.bam.bai",
            ref=config["genome"],
            ref_idx=get_fai(config["genome"]),
            ref_dict=get_dict(config["genome"]),
            known=config["dbsnp"],
            known_idx=get_tbi(config["dbsnp"])
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
            "/bio/gatk/baserecalibrator"

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

