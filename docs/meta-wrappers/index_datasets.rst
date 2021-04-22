.. _`index_datasets`:

INDEX_DATASETS
==============

Index genome sequences


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_index_datasets = {
        "genome": "reference/genome.fasta"
    }

    try:
        if config == dict():
            config = default_config_index_datasets
    except NameError:
        config = default_config_index_datasets



    def get_fasta_index_from_genome_path(genome_path: str) -> str:
        return genome_path + ".fai"


    def get_dict_from_genome_path(genome_path: str) -> str:
        return ".".join(genome_path.split(".")[:-1]) + ".dict"

    rule target:
        input:
            get_dict_from_genome_path(config["genome"]),
            get_fasta_index_from_genome_path(config["genome"])
        message:
            "Finishing indexing meta-wrapper"


    """
    This rule indexes the input genome sequence with Samtools. It is not
    explicitely requested by Samtools, but it will crash if the genome sequence
    is not indexed.

    This rule is cached since it should be used only once per reference sequence
    """
    rule samtools_faidx:
        input:
            config["genome"]
        output:
            get_fasta_index_from_genome_path(config["genome"])
        message: "Indexing reference fasta with Samtools"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 4098),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            ""
        log:
            "logs/samtools/faidx/genome.log"
        wrapper:
            "/bio/samtools/faidx"


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
            get_dict_from_genome_path(config["genome"])
        message: "Creating sequence dictionnary over reference genome with Picard"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 2048, 8192),
            time_min=lambda wildcards, attempt: attempt * 45
        params:
            ""
        log:
            "logs/picard/create_sequence_dictionnary/genome.log"
        wrapper:
            "/bio/picard/createsequencedictionary"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/samtools/faidx`

* :ref:`bio/picard/createsequencedictionary`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

