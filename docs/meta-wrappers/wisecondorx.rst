.. _`wisecondorx`:

WISECONDORX
===========

Search CNV in low pass WGS with WisecondorX


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_wisecondorx_config = {
        "convert_extra": "",
        "newref_extra": "",
        "predict_extra": "",
        "gender_extra": "",
        "ref": "/path/to/reference.fasta",
        "reference_samples": [],

    }

    rule wisecondorx_convert:
        input:
            aln = "bowtie2/map/{sample}.bam",
            ref = config["ref"]
        output:
            "wisecondorx/convert/{sample}.npz"
        message: "Converting reads {wildcards.sample} to NPZ"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            extra=config.get("convert_extra", "")
        log:
            "logs/wisecondorx/convert/{sample}.log"
        wrapper:
            "bio/wisecondorx/convert"


    rule wisecondorx_newref:
        input:
            aln = expand(
                "wisecondorx/convert/{sample}.npz",
                sample=config["reference_samples"]
            )
        output:
            "wisecondorx/reference.npz"
        message: "Building reference with WisecondorX"
        threads: 8
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            extra=config.get("newref_extra", "")
        log:
            "logs/wisecondorx/reference.log"
        wrapper:
            "bio/wisecondorx/newref"


    rule wisecondorx_predict:
        input:
            test = "wisecondorx/convert/{sample}.npz",
            ref = "wisecondorx/reference.npz"
        output:
            bins = "wisecondorx/{sample}/{sample}_bins.bed",
            segments = "wisecondorx/{sample}/{sample}_segments.bed",
            aberrations ="wisecondorx/{sample}/{sample}_aberrations.bed",
            statistics = "wisecondorx/{sample}/{sample}_statistics.bed"
        message: "Predict copy number alterations for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            extra = config.get("predict_extra", ""),
            prefix = lambda wildcards: f"wisecondorx/{wildcards.sample}/{wildcards.sample}"
        log:
            "logs/wisecondorx/predict/{sample}.log"
        wrapper:
            "bio/wisecondorx/predict"


    rule wisecondorx_gender:
        input:
            test = "wisecondorx/convert/{sample}.npz",
            ref = "wisecondorx/reference.npz"
        output:
            "wisecondorx/convert/{sample}.gender.txt"
        message: "Acessing gender of {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            extra = config.get("gender_extra", "")
        log:
            "logs/wisecondorx/gender/{sample}.log"
        wrapper:
            "bio/wisecondorx/gender"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/wisecondorx/convert`

* :ref:`bio/wisecondorx/newref`

* :ref:`bio/wisecondorx/predict`

* :ref:`bio/wisecondorx/gender`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

