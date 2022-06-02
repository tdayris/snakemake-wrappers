.. _`somatic_tmb`:

SOMATIC_TMB
===========

Compute tumor mutational burden within somatic samples


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config = {
        "filter_in": [],
        "filter_out": [],
        "min_coverage": 10,
        "allele_depth_keyname": "AD",
        "bed": "/path/to/regions.bed",
        "tmb_highness_threshold": 10
    }

    def get_best_igs():
        """Select cached IGS if available"""
        if config["bed"] == "/mnt/beegfs/database/bioinfo/Index_DB/captureKitDesigns/hg38/v5/sureselect_v5.padded.nochr.bed":
            return "igs.surselect_v5.yaml"
        return "igs.yaml"

    """
    Compute Itegrated Genome Size from the capture kit bed

    Several known beds have a cached igs
    """
    rule estimate_igs:
        input:
            bed = config["bed"]
        output:
            yaml = temp("igs.yaml")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/igs_estimation.log"
        wrapper:
            "bio/tmb/igs_estimation"


    rule estimate_igs_sureselect_v5:
        input:
            bed = "/mnt/beegfs/database/bioinfo/Index_DB/captureKitDesigns/hg38/v5/sureselect_v5.padded.nochr.bed"
        output:
            yaml = "igs.surselect_v5.yaml"
        threads: 1
        cache: True
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/igs_estimation.surselect_v5.log"
        wrapper:
            "bio/tmb/igs_estimation"

    """
    Compute statistics over somatic mutations' coverage
    """
    rule extract_somatic_mutations:
        input:
            vcf = "snpsift/fixed/{sample}.vcf.gz"
        output:
            yaml = temp("tmb/{sample}.yaml")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/extract_somatic_mutations/{sample}.log"
        params:
            filter_in = config.get("filter_in", []),
            filter_out = config.get("filter_out", []),
            min_coverage = config.get("min_coverage", 10),
            allele_depth = config.get("allele_depth_keyname", "AD")
        wrapper:
            "bio/tmb/extract_somatic"

    """
    Compute Tumor Molecular Burden from previous indexes and stats
    """
    rule compute_tmb:
        input:
            igs = get_best_igs(),
            samples = expand(
                "tmb/{sample}.yaml", sample=config["sample_list"]
            )
        output:
            tsv = "TMB.tsv"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/tmb.log"
        params:
            high_threshold = config.get("tmb_highness_threshold", 10)
        wrapper:
            "bio/tmb/compute_tmb"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/tmb/igs_estimation`

* :ref:`bio/tmb/extract_somatic`

* :ref:`bio/tmb/compute_tmb`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

Key names are parameters passed through the configurarion dictionnary




Authors
-------


* Thibault Dayris

