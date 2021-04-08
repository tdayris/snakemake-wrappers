.. _`snpeff_annotate`:

SNPEFF_ANNOTATE
===============

Annotate a raw VCF with SnpEff


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    """
    This rule annotated using previously donwloaded reference annotations.
    """
    rule snpeff:
        input:
            calls="calling/{sample}.vcf",
            db="snpeff/download/ebola_zaire"
        output:
            calls="snpeff/calls/{sample}.vcf",
            stats="snpeff/report/{sample}.html",
            csvstats="snpeff/csvstats/{sample}.csv"
        message: "Annotating {wildcards.sample} with SnpEff"
        threads: 3
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcard, attempt: attempt * 90
        params:
            extra=config.get("snpeff_extra", "")
        log:
            "logs/snpeff/annotate/{sample}.log"
        wrapper:
            "0.72.0-529-gcfac9685e/bio/snpeff/annotate"


    """
    This rule download SnpEff database, which will be used for further annotations.

    This rule is cached, since it should be done only once per reference genome.
    """
    rule snpeff_download:
        output:
            directory("snpeff/download/{reference}")
        message: "Downloading {wildcards.reference} database for SnpEff"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
            time_min=lambda wildcard, attempt: attempt * 120
        params:
            reference="{reference}"
        log:
            "logs/snpeff/download/{reference}.log"
        wrapper:
            "0.72.0-529-gcfac9685e/bio/snpeff/download"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/snpeff/download`

* :ref:`bio/snpeff/annotate`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

