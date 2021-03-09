.. _`DESeq2_post_process`:

DESEQ2_POST_PROCESS
===================

Plot various information based on DESeq2 results


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    import os.path

    from typing import Optional

    def rst_exists(path: str) -> Optional[str]:
        """
        If a provided rst path exists, then it is returned,
        otherwise "None" is returned.
        """
        return path if os.path.exists(path) else None

    volcanoplot_rst = rst_exists("reports/bioinfokit_volcanoplot.rst")
    bioinfokit_heatmap_rst = rst_exists("reports/bioinfokit_heatmap.rst")
    bioinfokit_maplot_rst = rst_exists("reports/bioinfokit_maplot.rst")


    """
    This rule builds a Volcanoplot from DESeq2 results
    """
    rule bioinfokit_volcanoplot:
        input:
            "deseq2/wald/Cond_compairing_B_vs_A.tsv"
        output:
            report(
                "bioinfokit/figures/volcanoplot.png",
                caption=volcanoplot_rst,
                category="Volcano-Plot"
            )
        message: "Plotting Volcanoplot with Bioinfokit"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 2048, 10240),
            time_min=lambda wildcards, attempt: attempt * 20
        params:
            read_csv={
                "header": 0,
                "index_col": None,
                #"usecols": ["Unnamed: 0", "log2FoldChange", "padj"]
            },
            volcano={
                "lfc": "log2FoldChange",
                "pv": "padj",
                "geneid": "Unnamed: 0",
                "lfc_thr": config.get("lfc_thr", (0.001, 0.001)),
                "pv_thr": config.get("pv_thr", (0.05, 0.05)),
                "gstyle": 2,
                "sign_line": True,
                "plotlegend": True
            }
        log:
            "logs/bioinfokit/volcanoplot.log"
        wrapper:
            "0.71.1-473-g5572d4839/bio/bioinfokit/volcanoplot"

    """
    This rule creates a sample clustered heatmap from the filtered-counts table
    """
    rule bioinfokit_sample_heatmap:
        input:
            "deseq2/filtered/filtered_counts.tsv"
        output:
            report(
                "bioinfokit/figures/sample_heatmap.png",
                caption=bioinfokit_heatmap_rst,
                category="Clustered Heatmap"
            )
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
            "0.71.1-473-g5572d4839/bio/bioinfokit/heatmap"


    """
    This rule creates a MA-plot from DESeq2 merged results
    """
    rule bioinfokit_maplot:
        input:
            "deseq2/filtered/merged.tsv"
        output:
            report(
                "bioinfokit/figures/maplot.png",
                caption=bioinfokit_maplot_rst,
                category="MA-Plot"
            )
        message: "Building MA-plot"
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
            "logs/bioinfokit/maplot.png"
        wrapper:
            "0.71.1-473-g5572d4839/bio/bioinfokit/maplot"


    """
    This rule merges and filters both DESeq2 counts and results for further graphs
    """
    rule filter_deseq2:
        input:
            wald_tsv = "deseq2/wald/Cond_compairing_B_vs_A.tsv",
            dst_tsv = "deseq2/dst/Cond_compairing_B_vs_A.tsv",
            gene2gene = "tximport/gene2gene.tsv"
        output:
            filtered_counts="deseq2/filtered/filtered_counts.tsv",
            filtered_deseq2="deseq2/filtered/filtered_deseq2.tsv",
            merged_table="deseq2/filtered/merged.tsv"
        message: "Filtering and merging DESeq2 results"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: attempt * 4096,
            time_min=lambda wildcard, attempt: attempt * 20
        log:
            "logs/deseq2/filter.log"
        wrapper:
            "0.71.1-473-g5572d4839/bio/pandas/deseq2_merge"



    """
    This rule build the conversion table from transcript to genes and their names.
    """
    rule gene_to_gene:
        input:
            gtf="refs/ensembl/chr21.gtf"
        output:
            gene2gene_large="tximport/gene2gene.tsv"
        message: "Building transcripts/genes conversion table"
        cache: True
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: attempt * 1536,
            time_min=lambda wildcard, attempt: attempt * 45
        params:
            gencode = True,
            header = True,
            positions = True
        log:
            "logs/tximport/tx2gene.log"
        wrapper:
            "0.71.1-473-g5572d4839/bio/gtf/tx2gene"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/bioinfokit/heatmap`

* :ref:`bio/bioinfokit/maplot`

* :ref:`bio/bioinfokit/volcanoplot`

* :ref:`bio/pandas/deseq2_merge`

* :ref:`bio/gtf/tx2gene`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

MultiQC report is based on configurations from in-house scripts.




Authors
-------


* Thibault Dayris

