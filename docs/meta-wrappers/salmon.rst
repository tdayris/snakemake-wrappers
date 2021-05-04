.. _`salmon`:

SALMON
======

Pseudo map and quantify your reads over transcritome with `Salmon <https://salmon.readthedocs.io/en/latest/>`_


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    default_salmon_config = {
        "genome": "/path/to/sequence.fasta",
        "transcriptome": "/path/to/transcriptome.fasta"
    }

    """
    This rule pseudo-map and quantifies your paired reads over the indexed
    reference.
    """
    rule salmon_quant_paired:
        input:
            r1="reads/{sample}_1.fq.gz",
            r2="reads/{sample}_2.fq.gz",
            index="salmon/index",
            gtf=config["gtf"]
        output:
            quant="salmon/pseudo_mapping/{sample}/quant.sf",
            lib="salmon/pseudo_mapping/{sample}/lib_format_counts.json",
            mapping=temp("salmon/bams/{sample}.bam")
        message: "Quantifying {wildcards.sample} with Salmon"
        threads: min(config.get("threads", 20), 20)
        resources:
            time_min=lambda wildcards, attempt: attempt * 60,
            mem_mb=lambda wildcards, attempt: (
                min(attempt * 5120 + 2048, 20480)
            )
        params:
            libType = config.get("salmon_libtype", "A"),
            extra = config.get("salmon_quant_extra", "--numBootstraps 100 --validateMappings --gcBias --seqBias --posBias")
        log:
            "logs/salmon/quant/{sample}.log"
        wrapper:
            "/bio/salmon/quant"


    """
    This rule shows how to inherit from paired quantification rule to quantify
    unpaired reads
    """
    use rule salmon_quant_paired as salmon_quant_unpaired with:
        input:
            r="reads/{sample}.fq.gz",
            index="salmon/index"


    """
    Index your transcriptome or gentrome file with Salmon in order to map your
    reads against this reference.

    This rule is cached since it should be used only once per reference genome.
    """
    rule salmon_index:
        input:
            sequences="salmon/decoy/gentrome.fasta",
            decoys="salmon/decoy/decoys.txt"
        output:
            index=directory("salmon/index")
        message: "Indexing transcriptome/gentrome sequences with Salmon"
        cache: True
        threads: min(config.get("threads", 20), 20)
        resources:
            time_min=lambda wildcards, attempt, input: (
                attempt * (120 if "decoys" in input.keys() else 45)
            ),
            mem_mb=lambda wildcards, attempt, input: (
                attempt * (15360 if "decoys" in input.keys() else 10240)
            )
        params:
            extra=config.get("salmon_index_extra", "--keepDuplicates --gencode")
        log:
            "logs/salmon/index.log"
        wrapper:
            "/bio/salmon/index"


    """
    This rule is optional in case you want to use decoy sequences within your
    transcriptome. See salmon documentation for more information.

    This rule is cached since it should be used only once per reference genome.
    """
    rule salmon_generate_decoy_sequence:
        input:
            transcriptome=config["transcriptome"],
            genome=config["genome"]
        output:
            decoys="salmon/decoy/decoys.txt",
            gentrome="salmon/decoy/gentrome.fasta"
        message: "Building gentrome and decoy sequences for Salmon"
        cache: True
        threads: 2
        resources:
            time_min=lambda wildcards, attempt: min(attempt * 20, 30),
            mem_mb=lambda wildcards, attempt: attempt * 512
        log:
            "logs/salmon/decoys.log"
        wrapper:
            "/bio/salmon/generate_decoy"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/salmon/generate_decoy`

* :ref:`bio/salmon/index`

* :ref:`bio/salmon/quant`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

