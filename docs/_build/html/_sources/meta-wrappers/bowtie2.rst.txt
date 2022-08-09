.. _`bowtie2`:

BOWTIE2
=======

Map your reads with Bowtie2


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_bowtie2 = {
        "ref": {"genome": "/path/to/reference.fasta"},
        "bowtie2_build_extra": "",
        "bowtie2_map_extra": ""
    }


    rule bowtie2_map:
        input:
            sample=expand(
                "reads/{sample}.{stream}.fq.gz",
                stream=["1", "2"],
                allow_missing=True
            ),
        output:
            temp("bowtie2/mapped/{sample}.bam")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 20,
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir="tmp"
        params:
            index="bowtie2/index/genome",
            extra=config.get("bowtie2_extra", "")
        log:
            "logs/bowtie2/map/{sample}.log"
        wrapper:
            "bio/bowtie2/align"


    rule bowtie2_build:
        input:
            reference=config["ref"]["genome"]
        output:
            temp(multiext(
                "bowtie2/index/genome",
                ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
            )),
        message: "Indexing {input} with Bowtie2"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 20,
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir="tmp"
        log:
            "logs/bowtie2_build/build.log"
        params:
            extra=config.get("bowtie2_build_extra", "")
        wrapper:
            "bio/bowtie2/build"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/bowtie2/build`

* :ref:`bio/bowtie2/align`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

