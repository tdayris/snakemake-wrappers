.. _`SnpEff_SnpSift (under development)`:

SNPEFF_SNPSIFT (UNDER DEVELOPMENT)
==================================

Annotate VCF files with SnpEff and SNpSift

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your directory containing VCF files

  cd /path/to/VCF/dir

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/snpeff_snpsift/run.sh hg38

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/snpeff_snpsift/run.sh hg19


Input/Output
------------


**Input:**

 
  
* VCF files (searched in current working directory)
  
 


**Output:**

 
  
   
* ``snpsift``: Annotated VCF file
   
  
 





Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:


* :ref:`meta/bio/snpeff_annotate`

* :ref:`meta/bio/snpsift`


Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.




Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`bio/tabix`

* :ref:`bio/compress/pbgzip`

* :ref:`bio/multiqc`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    import logging
    import os
    import pandas
    import sys
    from pathlib import Path

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../common/python")
    sys.path.append(common)

    from file_manager import *
    from files_linker import *
    from write_yaml import *
    from messages import *
    from snakemake.utils import min_version
    min_version("6.0")

    default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
    configfile: get_config(default_config)
    design = build_design(os.getcwd(), search_fastq_pairs)

    from pathlib import Path


    samples_list = design["Sample_id"]


    wildcard_constraints:
        sample = r"|".join(samples_list)


    last_vcf = (
        "snpsift/gnomad/{sample}.vcf"
        if config["params"]["NBCI_build"] != "mm10"
        else "snpsift/dbsnp/{sample}.vcf"
    )


    rule all:
        input:
            calls=expand(
                "snpsift/fixed/{sample}.vcf.gz",
                sample=samples_list,
                index=["", ".tbi"]
            ),
            qc="multiqc/SnpEff_annotation.html",
            tsv=expand(
                "snpsift/extractFields/{sample}.tsv",
                sample=samples_list
            )
        message:
            "Finishing the annotation pipeline"

    #################################
    ### FINAL VCF FILE INDEXATION ###
    #################################

    module compress_index_vcf_meta:
        snakefile: "../../meta/bio/compress_index_vcf/test/Snakefile"
        config: config


    use rule pbgzip_compress from compress_index_vcf_meta with:
        output:
            protected("{tool}/{subcommand}/{sample}.vcf.gz")


    use rule tabix_index from compress_index_vcf_meta with:
        output:
            protected("{tool}/{subcommand}/{sample}.vcf.gz.tbi")
        threads: 1


    #####################
    ### Export to TSV ###
    #####################

    rule extractfields:
        input:
            call="snpsift/fixed/{sample}.vcf.gz",
            call_index=get_tbi("snpsift/fixed/{sample}.vcf.gz")
        output:
            tsv=protected("snpsift/extractFields/{sample}.tsv")
        message:
            "Making {wildcards.sample} annotated VCF readable"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcards, attempt: attempt * 20
        log:
            "logs/snpsift/extractAllFields/{sample}.log"
        params:
            extra="-s '\t' -e '.'"
        wrapper:
            "bio/snpsift/extractAllFields"


    rule fix_vcf:
        input:
            vcf="snpsift/splitted/{sample}.vcf"
        output:
            vcf=temp("snpsift/fixed/{sample}.vcf")
        message:
            "Removing empty fields, trailing ';' and non-canonical chromosomes "
            "for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/bigr_scripts/fix_vcf/{sample}.log"
        params:
            default_chr=[*map(str, range(23)), *range(23), "MT", "X", "Y"],
            remove_non_conventional_chromosomes=True
        wrapper:
            "bio/BiGR/fix_vcf"


    rule split_vcf_features:
        input:
            call="snpsift/format2info/{sample}.vcf"
        output:
            call=temp("snpsift/splitted/{sample}.vcf")
        message:
            "Splitting ANN field of {wildcards.sample} to better "
            "end-user TSV experience"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/bigr_scripts/split_ann/{sample}.log"
        wrapper:
            "bio/BiGR/split_vcf_features"


    rule format_to_info:
        input:
            call="snpsift/gnomad/{sample}.vcf"
        output:
            call=temp("snpsift/format2info/{sample}.vcf")
        message:
            "Moving FORMAT data to INFO to increase readability and fixing field "
            "name redundancy for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/bigr_scripts/format2info/{sample}.log"
        wrapper:
            "bio/BiGR/vcf_format_to_info"


    ###############
    ### MultiQC ###
    ###############

    rule multiqc:
        input:
            expand(
                "snpeff/report/{sample}.html",
                sample=samples_list
            ),
            expand(
                "snpeff/csvstats/{sample}.csv",
                sample=samples_list
            )
        output:
            report(
                "multiqc/SnpEff_annotation.html",
                caption="../common/reports/multiqc.rst",
                category="Quality Controls"
            )
        message:
            "Aggregating quality reports from SnpEff"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        log:
            "logs/multiqc.log"
        wrapper:
            "bio/multiqc"


    ######################
    ### VCF annotation ###
    ######################

    snpeff_snpsift_config = {
        "ref": config["ref"],
        **config["snpeff_snpsift"]
    }

    module snpeff_meta:
        snakefile: "../../meta/bio/snpeff_annotate/test/Snakefile"
        config: snpeff_snpsift_config

    use rule snpeff from snpeff_meta with:
        input:
            calls="calls/{sample}.vcf.gz",
            calls_index="calls/{sample}.vcf.gz.tbi",
            db=config["ref"]["snpeff"]
        output:
            calls=temp("snpeff/{sample}.vcf"),
            stats=temp("snpeff/{sample}.html"),
            csvstats=temp("snpeff/{sample}.csv")


    module snpsift:
        snakefile: "../../meta/bio/snpsift/test/Snakefile"
        config: snpeff_snpsift_config

    use rule * from snpsift




Authors
-------


* Thibault Dayris
