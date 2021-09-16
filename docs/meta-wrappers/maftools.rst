.. _`MAFtools`:

MAFTOOLS
========

Plot various information about your samples and your cohort with MAFtools


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    default_config_maftools = {
        "genome": "GRCh38"
    }


    """
    Compute cosine similarity with known mutational signatures
    """
    rule maftools_cosine_similarity:
        pass

    """
    Compute signatures
    """
    rule maftools_signatures:
        pass

    """
    Build the trinucleotide matrix for signature analysis, on GRCh27
    """
    rule maftools_trinucleotidematrix_hg19:
        input:
            rds="maf/maftools_samples/{sample}/maf.RDS"
        output:
            tsv="maf/maftools_samples/{sample}/trinucleotide_matrix_GRCh27.tsv",
            #png="maf/maftools_samples/{sample}/trinucleotide_signtures_hg19.tsv",
            rds="maf/maftools_samples/{sample}/trinucleotide_matrix_GRCh27.RDS",
        message: "Computing trinucleotide matrices for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/maftools/trinucleotide_matrix/{sample}.log"
        params:
            estimate_extra=config.get("estimate_extra", "nTry=3")
        wrapper:
            "bio/maftools/trinucleotidematrix_hg19"


    """
    Build the trinucleotide matrix for signature analysis, on GRCh38
    """
    rule maftools_trinucleotidematrix_hg38:
        input:
            rds="maf/maftools_samples/{sample}/maf.RDS"
        output:
            tsv="maf/maftools_samples/{sample}/trinucleotide_matrix_hg38.tsv",
            #png="maf/maftools_samples/{sample}/trinucleotide_signtures_hg38.tsv",
            rds="maf/maftools_samples/{sample}/trinucleotide_matrix_hg38.RDS",
        message: "Computing trinucleotide matrices for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/maftools/trinucleotide_matrix/{sample}.log"
        params:
            estimate_extra=config.get("estimate_extra", "nTry=3")
        wrapper:
            "bio/maftools/trinucleotidematrix_hg38"


    """
    Load a maf file into R for further MAFtools analysis
    """
    rule maftools_readmaf:
        input:
            "maf/translated/{sample}.maf"
        output:
            rds="maf/maftools_samples/{sample}/maf.RDS",
            summary=temp("maf/maftools_samples/{sample}/maf_summary.txt")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            summay_prefix=lambda w: f"maf/maftools_samples/{w.sample}/maf"
        wrapper:
            "bio/maftools/readmaf"


    """
    SnpEff contains multiple variants effects, given the wide range of biology
    possibilities. MAFtools keeps only a set of six classes, or its graphs would
    become unreadable.
    This rule makes the translation
    """
    rule maftools_rename_snpeff_effect:
        input:
            maf="maf/maftools/{sample}.maf"
        output:
            maf=temp("maf/translated/{sample}.maf")
        message: "Changing snpeff effects in {wildcards.sample} for MAFtools"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/maftools/translated/{sample}.log"
        wrapper:
            "bio/maftools/translate_snpeff_maf_effect"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/maftools/translate_snpeff_maf_effect`

* :ref:`bio/maftools/readmaf`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

No inclusion of clinical data possible yet.




Authors
-------


* Thibault Dayris

