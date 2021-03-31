.. _`bwa_fixmate`:

BWA_FIXMATE
===========

Map your reads with BWA and fix mates right after


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    try:
        if config == dict():
            config = {"threads": 20}
    except NameError:
        config = {"threads": 20}

    """
    This rule indexes the bam file since almost all downstream tools requires it
    """
    rule samtools_index:
        input:
            "samtools/sort/{sample}.bam"
        output:
            "samtools/sort/{sample}.bam.bai"
        message: "Indexing mapped reads of {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=1536,
            time_min=lambda wildcards, attempt: attempt * 45
        log:
            "logs/samtools/sort/{sample}.log"
        wrapper:
            "0.72.0-504-g15bd6157e/bio/samtools/index"


    """
    This rule sorts reads by position for further analyses.
    """
    rule samtools_sort_coordinate:
        input:
            "samtools/fixmate/{sample}.bam"
        output:
            "samtools/sort/{sample}.bam"
        message:
            "Sorting {wildcards.sample} reads by query name for fixing mates"
        threads: config.get("threads", 20)
        resources:
            mem_mb=lambda wildcards, threads: threads * 1792,
            time_min=lambda wildcards, attempt: attempt * 90
        log:
            "logs/samtools/query_sort_{sample}.log"
        params:
            extra = "-m 1536M"
        wrapper:
            "0.72.0-504-g15bd6157e/bio/samtools/sort"


    """
    BWA sometimes fails to annotate read mates correctly. We fix this behaviour
    with the rule below.
    """
    rule samtools_fixmate:
        input:
            "bwa_mem2/mem/{sample}.bam"
        output:
            temp("samtools/fixmate/{sample}.bam")
        message: "Fixing mate annotation on {wildcards.sample} with Samtools"
        threads: config.get("threads", 20)
        resources:
            mem_mb = (
                lambda wildcards, attempt: min(attempt * 2048 + 2048, 8192)
            ),
            time_min = (
                lambda wildcards, attempt: min(attempt * 45, 180)
            )
        params:
            config.get("fixmate_extra", "-cmr")
        log:
            "logs/samtools/fixmate/{sample}.log"
        wrapper:
            "0.72.0-504-g15bd6157e/bio/samtools/fixmate"


    """
    This rule maps your reads against the indexed reference with BWA.
    """
    rule bwa_mem:
        input:
            reads = expand(
                "reads/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
            index=multiext(
                "bwa_mem2/index/genome", ".0123", ".amb", ".ann", ".pac"
            )
        output:
            temp("bwa_mem2/mem/{sample}.bam")
        message: "Mapping {wildcards.sample} with BWA"
        threads: config.get("threads", 20)
        resources:
            mem_mb = (
                lambda wildcards, attempt: min(attempt * 6144 + 2048, 20480)
            ),
            time_min = (
                lambda wildcards, attempt: min(attempt * 120, 480)
            )
        params:
            index="bwa_mem2/index/genome",
            extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
            sort="samtools",         # We chose Samtools to sort by queryname
            sort_order="queryname",  # Queryname sort is needed for a fixmate
            sort_extra="-m 1536M"     # We extand the sort buffer memory
        log:
            "log/bwa_mem2/mem/{sample}.log"
        wrapper:
            "0.72.0-504-g15bd6157e/bio/bwa-mem2/mem"


    """
    Index your reference genome with BWA.

    This rule is cached since it should be used once per reference genome
    """
    rule bwa_index:
        input:
            "sequence/genome.fasta"
        output:
            multiext(
                "bwa_mem2/index/genome", ".0123", ".amb", ".ann", ".pac"
            )
        message: "Indexing reference genome with BWA"
        cache: True
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: min(attempt * 90, 480),
            mem_mb=lambda wildcards, attempt: min(attempt * 6144 + 2048, 20480)
        params:
            prefix="bwa_mem2/index/genome"
        log:
            "logs/bwa_mem2/index/genome.log"
        wrapper:
            "0.72.0-504-g15bd6157e/bio/bwa-mem2/index"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/samtools/index`

* :ref:`bio/samtools/sort`

* :ref:`bio/samtools/fixmate`

* :ref:`bio/bwa-mem2/mem`

* :ref:`bio/bwa-mem2/index`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

