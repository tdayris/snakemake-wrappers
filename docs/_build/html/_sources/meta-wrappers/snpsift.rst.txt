.. _`snpsift`:

SNPSIFT
=======

Annotate SNP/Indels with snpsift


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config={
        "samples":["test"],
        "ref": {
            "cosmic": "/path/to/annotation",
            "dbsnp": "/path/to/annotation",
            "dbnsfp": "/path/to/annotaion",
            "fasta": "/path/to/annotation",
            "gmt": "/path/to/annotation",
            "gwascat": "/path/to/annotation",
            "kaviar": "/path/to/annotation",
        }
    }

    try:
        if config == dict():
            config = default_config
    except NameError:
        config = default_config


    #rule all:
    #    input:
    #        expand(
    #            "snpsift/cosmic/{sample}.vcf",
    #            sample=config["samples"]
    #        )


    rule snpsift_dbnsfp:
        input:
            call = "snpsift/gwascat/{sample}.vcf",
            dbNSFP = config["ref"]["dbnsfp"],
            dbNSFP_tbi = config["ref"]["dbnsfp"] + ".tbi"
        output:
            call = temp("snpsift/dbnsfp/{sample}.vcf")
        message:
            "Annotating {wildcards.sample} with dbNSFP"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/snpsift/dbnsfp/{sample}.log"
        wrapper:
            "/bio/snpsift/dbnsfp"



    rule snpsift_gwascat:
        input:
            call = "snpsift/cosmic/{sample}.vcf",
            gwascat = config["ref"]["gwascat"]
        output:
            call = temp("snpsift/gwascat/{sample}.vcf")
        message:
            "Annotating {wildcards.sample} with GWAS Catalog"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/snpsift/gwascat/{sample}.log"
        wrapper:
            "/bio/snpsift/gwascat"


    rule snpsift_cosmic:
        input:
            call="snpsift/dbsnp/{sample}.vcf",
            database=config["ref"]["cosmic"]
        output:
            call=temp("snpsift/cosmic/{sample}.vcf")
        message:
            "Annotating {wildcards.sample} with COSMIC"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/snpsift/cosmic/{sample}.log"
        wrapper:
            "/bio/snpsift/annotate"


    rule snpsift_dbsnp:
        input:
            call="snpsift/kaviar/{sample}.vcf",
            database=config["ref"]["dbsnp"]
        output:
            call=temp("snpsift/dbsnp/{sample}.vcf")

        message:
            "Annotating {wildcards.sample} with dbSNP"
        threads: 1
        log:
            "logs/snpsift/dbsnp/{sample}.log"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
            time_min=lambda wildcards, attempt: attempt * 45
        wrapper:
            "/bio/snpsift/annotate"


    rule snpsift_kaviar:
        input:
            call="snpsift/gmt/{sample}.vcf",
            database=config["ref"]["kaviar"]
        output:
            call=temp("snpsift/kaviar/{sample}.vcf")

        message:
            "Annotating {wildcards.sample} with Kaviar"
        threads: 1
        log:
            "logs/snpsift/kaviar/{sample}.log"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
            time_min=lambda wildcards, attempt: attempt * 45
        wrapper:
            "/bio/snpsift/annotate"


    rule snpsift_gmt:
        input:
            call = "snpsift/vartype/{sample}.vcf",
            gmt = config["ref"]["gmt"]
        output:
            call = temp("snpsift/gmt/{sample}.vcf")
        message:
            "Annotating {wildcards.sample} with MSigDB"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
            time_min=lambda wildcards, attempt: attempt * 45
        wrapper:
            "/bio/snpsift/genesets"


    rule snpsift_vartype:
        input:
            vcf="snpeff/calls/{sample}.vcf.gz",
            vcf_tbi="snpeff/calls/{sample}.vcf.gz.tbi"
        output:
            vcf=temp("snpsift/vartype/{sample}.vcf")
        message:
            "Annotating variant types in {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020 + 4096,
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/snpsift/varType/{sample}.log"
        wrapper:
            "/bio/snpsift/varType"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/snpsift/varType`

* :ref:`bio/snpsift/genesets`

* :ref:`bio/snpsift/annotate`

* :ref:`bio/snpsift/gwascat`

* :ref:`bio/snpsift/dbnsfp`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

