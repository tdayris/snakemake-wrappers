.. _`compress_index_vcf`:

COMPRESS_INDEX_VCF
==================

Compress and index VCF files


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    """
    This rule indexes the pbgzipped vcf file
    """
    rule tabix_index:
        input:
            "{tool}/{subcommand}/{sample}.vcf.gz"
        output:
            "{tool}/{subcommand}/{sample}.vcf.gz.tbi"
        message:
            "Indexing {wildcards.sample} calls with Tabix "
            "(built with {wildcards.tool} {wildcards.subcommand})"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 10240),
            time_min=lambda wildcards, attempt: attempt * 60
        params:
            "-p vcf"
        log:
            "logs/{tool}/{subcommand}/tabix/index/{sample}.log"
        wrapper:
            "bio/tabix"


    """
    This rule compress a provided VCF file with pbgzip
    """
    rule pbgzip_compress:
        input:
            "{tool}/{subcommand}/{sample}.vcf"
        output:
            "{tool}/{subcommand}/{sample}.vcf.gz"
        message:
            "Compressing {wildcards.sample} VCF file, "
            "(built with {wildcards.tool} {wildcards.subcommand})"
        threads:
            config.get("threads", 20)
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 10240),
            time_min=lambda wildcards, attempt: attempt * 30
        params:
            ""
        log:
            "logs/{tool}/{subcommand}/pbgzip/{sample}.log"
        wrapper:
            "bio/compress/pbgzip"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/tabix`

* :ref:`bio/compress/pbgzip`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

The VCF file must be in directories, which tree is formatted as follows: {tool}/{subcommand}/{sample}.vcf




Authors
-------


* Thibault Dayris

