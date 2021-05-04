.. _`meta_caller_germline`:

META_CALLER_GERMLINE
====================

Gather calls from Mutect2 and Varscan2


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    def get_fai(genome_path: str) -> str:
        return genome_path + ".fai"

    def get_tbi(vcf_path: str) -> str:
        return vcf_path + ".tbi"

    rule bcftools_merge:
        input:
            calls=[
                "varscan2/mpileup2cns_reheaded/{sample}.vcf.gz",
                "mutect2/filter_reheaded/{sample}.vcf.gz"
            ],
            calls_index=[
                get_tbi("varscan2/mpileup2cns_reheaded/{sample}.vcf.gz"),
                get_tbi("mutect2/filter_reheaded/{sample}.vcf.gz")
            ],
            regions=config["bed"]
        output:
            temp("meta_caller/calls/{sample}.vcf.gz")
        message:
            "Merging calls from Mutect2 and Varscan2 for {wildcards.sample}"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 8192, 20480),
            time_min=lambda wildcards, attempt: attempt * 30
        log:
            "logs/bcftools/merge/{sample}.log"
        params:
            "--merge none"
        wrapper:
            "/bio/bcftools/merge"

    ########################
    ### Renaming samples ###
    ########################

    rule rehead_before_merge:
        input:
            vcf="{tool}/{subcommand}/{sample}.vcf.gz",
            vcf_index=get_tbi("{tool}/{subcommand}/{sample}.vcf.gz"),
            fai=get_fai(config["genome"]),
            samples = "{tool}/{subcommand}/reheader/{sample}.renamed.sample.list",
        output:
            temp("{tool}/{subcommand}_reheaded/{sample}.vcf.gz")
        message:
            "Adding sequence info in header {wildcards.sample} "
            "(considering {wildcards.tool}, {wildcards.subcommand})"
        group:
            "Rename_{wildcards.sample}"
        threads: 2
        resources:
            mem_mb=512,
            time_min=10
        log:
            "logs/{tool}/{subcommand}/rehead/{sample}.log"
        wrapper:
            "/bio/bcftools/reheader"


    rule bcftools_name_correspondancy:
        input:
            "{tool}/{subcommand}/reheader/{sample}.sample.list",
            "{tool}/{subcommand}/reheader/caller_name.lst"
        output:
            temp("{tool}/{subcommand}/reheader/{sample}.renamed.sample.list")
        message:
            "Renaming sample list from {wildcards.sample}"
        group:
            "Rename_{wildcards.sample}"
        threads: 1
        resources:
            mem_mb=128,
            time_min=2
        log:
            "logs/{tool}/{subcommand}/rehead/{sample}.sample_list.renamed.log"
        shell:
            "paste {input} > {output} 2> {log}"


    rule bcftools_sample_list:
        input:
            "{tool}/{subcommand}/{sample}.vcf.gz"
        output:
            sample_list=temp("{tool}/{subcommand}/reheader/{sample}.sample.list")
        message:
            "Extracting sample list from {wildcards.sample}"
        group:
            "Rename_{wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 512, 4096),
            time_min=lambda wildcards, attempt: attempt * 5
        log:
            "logs/{tool}/{subcommand}/rehead/{sample}.sample_list.log"
        wrapper:
            "/bio/bcftools/query"


    rule tool_name:
        output:
            temp("{tool}/{subcommand}/reheader/caller_name.lst")
        message:
            "Building list of caller to rehead final vcf files"
        threads: 1
        resources:
            mem_mb=125,
            time_min=2
        log:
            "logs/{tool}/{subcommand}/rehead/tool_name.log"
        params:
            lambda wildcards: {
                "varscan2_mpileup2cns": "Varscan2",
                "mutect2_filter": "Mutect2",
            }.get(f"{wildcards.tool}_{wildcards.subcommand}")
        shell:
            "echo {params} > {output} 2> {log}"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/bcftools/merge`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

Input VCF files are supposed to be filtered and not annotated (merging annotations may lead to errors with BCFTools).




Authors
-------


* Thibault Dayris

