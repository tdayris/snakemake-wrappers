.. _`index_datasets`:

INDEX_DATASETS
==============

Index genome sequences


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    import sys
    from pathlib import Path

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../../../../bigr_pipelines/common/python")
    sys.path.append(common)

    from file_manager import *

    default_config_index_datasets = {
        "genome": "reference/genome.fasta"
    }

    try:
        if config == dict():
            config = default_config_index_datasets
    except NameError:
        config = default_config_index_datasets


    rule target:
        input:
            get_dict(config["genome"]),
            get_fai(config["genome"])
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
            get_fai(config["genome"])
        message: "Indexing reference fasta with Samtools"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 4098),
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir=lambda wildcards, input: f"tmp/genome_indexation.tmp"
        params:
            ""
        log:
            "logs/samtools/faidx/genome.log"
        wrapper:
            "bio/samtools/faidx"


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
            get_dict(config["genome"])
        message: "Creating sequence dictionnary over reference genome with Picard"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 2048, 8192),
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir=lambda wildcards, input: f"tmp/genome_indexation.tmp"
        params:
            ""
        log:
            "logs/picard/create_sequence_dictionnary/genome.log"
        wrapper:
            "bio/picard/createsequencedictionary"


    """
    This rule indexes a given bam file, using samtools index. Most of the time,
    fasta indexes are not explicitely requested by softwares, but they will crash
    if this index is missing.
    """
    rule samtools_index_bam:
        input:
            "{bam_path}.bam"
        output:
            "{bam_path}.bam.bai"
        message:
            "Indexing {bam_path} with Samtools index"
        threads: 20
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1024, 8192),
            time_min=lambda wildcards, attempt: attempt * 30,
            tmpdir=lambda wildcards: f"tmp/samtools_index.tmp"
        log:
            "samtools/{bam_path}/index.log"
        wrapper:
            "bio/samtools/index"

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

