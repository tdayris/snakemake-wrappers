.. _`meta_caller`:

META_CALLER
===========

Use multiple callers to assess variants in a single sample


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    from pathlib import Path
    from snakemake.utils import min_version
    min_version("6.0")


    def search_bam_samples(path: str):
        for f in Path(path).iterdir():
            if f.name.endswith(".bam") and f.is_file():
                yield(str(f.name)[:-len(".bam")])
            elif f.is_dir() and not f.name.startswith(".snakemake"):
                yield from search_bam_samples(f.absolute())



    def get_fasta_index_from_genome_path(genome_path: str) -> str:
        return genome_path + ".fai"


    def get_dict_from_genome_path(genome_path: str) -> str:
        return ".".join(genome_path.split(".")[:-1]) + ".dict"


    default_config_meta_caller = {
        "samples": list(search_bam_samples(".")),
        "genome": "reference/genome.fasta",
        "threads": 20
    }


    try:
        if config == dict():
            config = default_config_meta_caller
    except NameError:
        config = default_config_meta_caller


    tools_and_subcommands = [
        ("varscan2", "concat"),
        ("gatk", "mutect2"),
        ("strelka", "germline")
    ]

    wildcard_constraints:
        sample = r"|".join(config["samples"]),
        tool = r"|".join([i[0] for i in tools_and_subcommands]),
        subcommand = r"|".join([i[-1] for i in tools_and_subcommands])


    ruleorder: index_datasets_samtools_faidx > varscan2_calling_samtools_faidx


    rule target:
        input:
            expand("bcftools/merge/{sample}.vcf.gz", sample=config["samples"])
        message:
            "Finishing the meta_caller meta-wrapper"


    ################################
    ###     Modules loading      ###
    ################################


    module varscan2_calling:
        snakefile: "../../varscan2_calling/test/Snakefile"
        config: config

    module index_datasets:
        snakefile: "../../index_datasets/test/Snakefile"
        config: config


    use rule * from varscan2_calling as varscan2_calling_*

    use rule * from index_datasets as index_datasets_*

    ################################
    ### Merge all separate calls ###
    ################################


    rule bcftools_merge:
        input:
            calls=[
                "varscan2/concat/{sample}.fai.vcf.gz",
                "gatk/mutect2/{sample}.fai.vcf.gz",
                "strelka/germline/{sample}.fai.vcf.gz"
            ],
            calls_indexes=[
                "varscan2/concat/{sample}.fai.vcf.gz.tbi",
                "gatk/mutect2/{sample}.fai.vcf.gz.tbi",
                "strelka/germline/{sample}.fai.vcf.gz.tbi"
            ],
            strelka_results = "strelka/germline/{sample}"
        output:
            "bcftools/merge/{sample}.vcf.gz"
        message:
            "Merging multiple calls for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/bcftools/merge/{sample}.log"
        params:
            "-m none --force-samples --missing-to-ref --no-index"
        wrapper:
            "/bio/bcftools/merge"


    rule rehead_before_merge:
        input:
            vcf="{tool}/{subcommand}/{sample}.vcf.gz",
            fai=get_fasta_index_from_genome_path(config["genome"]),
            samples = "{tool}/{subcommand}/reheader/{sample}.renamed.sample.list",
        output:
            temp("{tool}/{subcommand}/{sample}.fai.vcf.gz")
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
                "varscan2_concat": "Varscan2",
                "gatk_mutect2": "Mutect2",
                "strelka_germline": "Strelka"
            }.get(f"{wildcards.tool}_{wildcards.subcommand}")
        shell:
            "echo {params} > {output} 2> {log}"


    rule pbgzip:
        input:
            "{tool}/{subcommand}/{sample}.vcf"
        output:
            temp("{tool}/{subcommand}/{sample}.vcf.gz")
        message:
            "Compressing {sample} calling with Prallel Block GZip "
            "(for {wildcards.tool}, {wildcards.subcommand})"
        threads: config.get("threads", 20)
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 512, 2048),
            time_min=lambda wildcards, attempt: attempt * 30
        log:
            "logs/{tool}/{subcommand}/pbgzip/{sample}.log"
        wrapper:
            "/bio/compress/pbgzip"


    ########################
    ### Varscan2 Calling ###
    ########################


    use rule bcftools_concat from varscan2_calling with:
        output:
            "varscan2/concat/{sample}.vcf.gz"


    use rule tabix_index from varscan2_calling with:
        input:
            "{tool}/{subcommand}/{sample}.vcf.gz"
        output:
            temp("{tool}/{subcommand}/{sample}.vcf.gz.tbi")
        message:
            "Indexing compressed VCF for sample {wildcards.sample} "
            "(for {wildcards.tool}, {wildcards.subcommand})"
        log:
            "logs/{tool}/{subcommand}/tabix/index/{sample}.log"


    use rule varscan2_calling_tabix_index as tabix_index_post_fai with:
        input:
            "{tool}/{subcommand}/{sample}.fai.vcf.gz"
        output:
            "{tool}/{subcommand}/{sample}.fai.vcf.gz.tbi"
        log:
            "logs/{tool}/{subcommand}/tabix/index/{sample}.fai.log"


    ####################
    ### GATK Calling ###
    ####################


    rule mutect2:
        input:
            fasta = config["genome"],
            fasta_index = get_fasta_index_from_genome_path(config["genome"]),
            fasta_dict = get_dict_from_genome_path(config['genome']),
            map = "mapped/{sample}.bam",
            map_index = "mapped/{sample}.bam.bai",
        output:
            vcf = "gatk/mutect2/{sample}.vcf.gz"
        message:
            "Calling variants in {wildcards.sample} with Mutect2"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 8192, 20480),
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/gatk/mutect/{sample}.log"
        wrapper:
             "/bio/gatk/mutect"


    #######################
    ### Strelka Calling ###
    #######################


    rule strelka_rename:
        input:
            "strelka/germline/{sample}"
        output:
            vcf="strelka/germline/{sample}.vcf.gz",
            tbi="strelka/germline/{sample}.vcf.gz.tbi"
        message:
            "Renaming strelka output"
        threads: 1
        resources:
            mem_mb=128,
            time_min=5
        log:
            "logs/strelka/rename/{sample}.log"
        params:
            "-sfr"
        shell:
            "ln {params} {input}/results/variants/variants.vcf.gz {output.vcf} "
            "> {log} 2>&1 && "
            "ln {params} {input}/results/variants/variants.vcf.gz.tbi {output.tbi} "
            ">> {log} 2>&1 "


    rule strelka_germline:
        input:
            bam = "mapped/{sample}.bam",
            bam_index = "mapped/{sample}.bam.bai",
            fasta = config["genome"],
            fasta_index = get_fasta_index_from_genome_path(config["genome"])
        output:
            temp(directory("strelka/germline/{sample}"))
        message:
            "Running strelka on {wildcards.sample}"
        threads: config.get("threads", 20)
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 8192, 2048),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            run_extra = "",
            config_extra = ""
        log:
            "logs/strelka/{sample}.log"
        wrapper:
            "/bio/strelka/germline"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/bcftools/concat`

* :ref:`bio/samtools/mpileup`

* :ref:`bio/tabix`

* :ref:`bio/varscan/mpileup2indel`

* :ref:`bio/varscan/mpileup2snp`

* :ref:`bio/gatk/mutect2`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

This meta caller uses:

* Varscan2 (snp + indel)
* Mutect2




Authors
-------


* Thibault Dayris

