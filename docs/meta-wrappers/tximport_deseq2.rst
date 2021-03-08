.. _`tximport_deseq2`:

TXIMPORT_DESEQ2
===============

Import counts with tximport and run DESeq2


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    """
    This rule performs the size factor and dispersions estimations as well as the
    wald test.
    """
    rule deseq2:
        input:
            dds="deseq2/dds.RDS"
        output:
            rds = "deseq2/wald/Cond_compairing_B_vs_A.RDS",
            deseq2_result_dir = directory("deseq2/tsv"),
            deseq2_tsv = "deseq2/wald/Cond_compairing_B_vs_A.tsv",
            normalized_counts = "deseq2/dst/counts.tsv",
            dst = "deseq2/dst/counts.RDS"
        message: "Running DESeq2 analysis"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: attempt * 2048,
            time_min=lambda wildcard, attempt: attempt * 60
        params:
            contrast = config.get("deseq2_contrast", ["Cond", "B", "A"]),
        log:
            "logs/deseq2/deseq.log"
        wrapper:
            "0.71.1-460-g2d0d5bf6e/bio/deseq2/DESeq"


    """
    This rule formats counts for DESeq2. The design matrix and its corresponding
    formula are included.
    """
    rule deseq2_dataset_from_tximport:
        input:
            tximport="tximport/txi.RDS",
            coldata="coldata.tsv",
        output:
            dds="deseq2/dds.RDS"
        message: "Formatting counts for DESeq2",
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 3072, 20480),
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            design=config.get("deseq2_formula", "~Condition"),
            levels=config.get("deseq2_levels", ["A", "B"]),
            factor=config.get("deseq2_factor", "Condition")
        log:
            "logs/deseq2/deseq2_dataset_from_tximport.log"
        wrapper:
            "0.71.1-460-g2d0d5bf6e/bio/deseq2/DESeqDataSetFromTximport/"


    """
    This rule imports counts from tables to R data object. Its memory requirements
    are linked to the number of samples
    """
    rule tximport:
        input:
            quant=expand(
                "quant/{sample}/quant.sf",
                sample=[f"{chr(i)}.chr21" for i in range(97, 103)]
            )
        output:
            txi=temp("tximport/txi.RDS")
        message: "Importing counts in DESeq2"
        threads: 1
        resources:
            mem_mb=lambda wildcard, input: len(input.quant) * 1024,
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            extra=config.get("tximport_extra", "type='salmon'")
        log:
            "logs/tximport.log"
        wrapper:
            "0.71.1-460-g2d0d5bf6e/bio/tximport"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/tximport`

* :ref:`bio/deseq2/DESeqDataSetFromTximport`

* :ref:`bio/deseq2/DESeq`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

The R statistical formula must refer to columns in design file.




Authors
-------


* Thibault Dayris

