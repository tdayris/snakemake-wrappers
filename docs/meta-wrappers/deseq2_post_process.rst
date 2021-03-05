.. _`DESeq2_post_process`:

DESEQ2_POST_PROCESS
===================

Plot various information based on DESeq2 results


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule bioinfokit_volcanoplot:
        input:
            "deseq2/wald/results.tsv"
        output:
            "bioinfokit/volcanoplot.png"
        message: "Plotting Volcanoplot with Bioinfokit"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 2048, 10240),
            time_min=lambda wildcards, attempt: attempt * 20
        params:
            read_csv={
                "header": 0,
                "index_col": None,
                "usecols": ["Unnamed: 0", "log2FoldChange", "padj"]
            }
            volcano={
                "lfc": "log2FoldChange",
                "pv": "padj",
                "geneid": "Unnamed: 0",
                "gstyle": 2,
                "sign_line": True
            }
        log:
            "logs/bioinfokit/volcanoplot.log"
        wrapper:
            "0.71.1-459-g6aed01be9/bio/bioinfokit/volcanoplot"


    rule bioinfokit_sample_heatmap:
        input:
            "deseq2/dst/counts.tsv"
        output:
            "bioinfokit/sample_heatmap.png"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 2048, 10240),
            time_min=lambda wildcards, attempt: attempt * 20
        params:
            read_csv={
                "header": 0,
                "index_col": 0
            },
            hmap={
                "rowclus": True,
                "colclus": True
            }
        log:
            "logs/bioinfokit/sample_heatmap.png"
        wrapper:
            "0.71.1-459-g6aed01be9/bio/bioinfokit/heatmap"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/bioinfokit/volcanoplot`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

MultiQC report is based on configurations from in-house scripts.




Authors
-------


* Thibault Dayris

