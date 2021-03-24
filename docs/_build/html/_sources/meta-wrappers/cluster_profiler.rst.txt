.. _`cluster_profiler`:

CLUSTER_PROFILER
================

Perform gene sets analysis against pathways


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule cluster_profiler_upsetplot:
        input:
            genelist="clusterProfiler/genelist.RDS"
        output:
            png="clusterProfiler/pathview/{pathway_id}.png"
        message: "Building pathview over {wildcards.pathway_id}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        params:
            pathway_id="{pathway_id}"
        log:
            "logs/clusterProfiler/pathview/{pathway_id}.log"
        wrapper:
            "0.72.0-493-g8b815973b/bio/clusterProfiler/pathview"


    rule cluster_profiler_upsetplot:
        input:
            rds="clusterProfiler/{database}/enrich.RDS"
        output:
            png="clusterProfiler/{database}/upsetplot.png"
        message: "Building upsetplot over {wildcards.database}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        params:
            upsetplot_extra=config.get("upsetplot_extra", "")
        log:
            "logs/clusterProfiler/upsetplot/upsetplot_{database}.log"
        wrapper:
            "0.72.0-493-g8b815973b/bio/clusterProfiler/upsetplot"


    rule cluster_profiler_heatplot:
        input:
            rds="clusterProfiler/{database}/enrich.RDS",
            gene_list="clusterProfiler/genelist.RDS"
        output:
            png="clusterProfiler/{database}/heatplot.png"
        message: "Building heatplot over {wildcards.database}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        params:
            heatplot_extra=config.get("heatplot_extra", "")
        log:
            "logs/clusterProfiler/heatplot/heatplot_{database}.log"
        wrapper:
            "0.72.0-493-g8b815973b/bio/clusterProfiler/heatplot"


    rule cluster_profiler_cnetplot:
        input:
            rds="clusterProfiler/{database}/enrich.RDS",
            gene_list="clusterProfiler/genelist.RDS"
        output:
            png="clusterProfiler/{database}/cnetplot.png"
        message: "Building cnetplot over {wildcards.database}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        params:
            cnetplot_extra=config.get("cnetplot_extra", "")
        log:
            "logs/clusterProfiler/cnetplot/cnetplot_{database}.log"
        wrapper:
            "0.72.0-493-g8b815973b/bio/clusterProfiler/cnetplot"


    rule cluster_profiler_dotplot:
        input:
            rds="clusterProfiler/{database}/enrich.RDS"
        output:
            png="clusterProfiler/{database}/dotplot.png"
        message: "Building dotplot over {wildcards.database}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        params:
            dotplot_extra=config.get("dotplot_extra", "")
        log:
            "logs/clusterProfiler/dotplot/dotplot_{database}.log"
        wrapper:
            "0.72.0-493-g8b815973b/bio/clusterProfiler/dotplot"


    rule cluster_profiler_barplot:
        input:
            rds="clusterProfiler/{database}/enrich.RDS"
        output:
            png="clusterProfiler/{database}/barplot.png"
        message: "Building barplot over {wildcards.database}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        params:
            barplot_extra=config.get("barplot_extra", "")
        log:
            "logs/clusterProfiler/barplot/barplot_{database}.log"
        wrapper:
            "0.72.0-493-g8b815973b/bio/clusterProfiler/barplot"


    rule cluster_profiler_enrich_go:
        input:
            rds="clusterProfiler/genelist.RDS"
        output:
            rds="clusterProfiler/GO_{onthology}/enrich.RDS",
            tsv="clusterProfiler/GO_{onthology}/enrichGO_{onthology}.tsv"
        message: "Running GO:{wildcards.onthology} enrichment"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        params:
            extra=lambda wildcards: f" ont={wildcards.onthology} "
        log:
            "logs/clusterProfiler/enrich/GO_{onthology}/enrichGO_{onthology}.log"
        wrapper:
            "0.72.0-493-g8b815973b/bio/clusterProfiler/enrichGO"


    rule translate_genes_identifiers_hg38:
        input:
            table="deseq2/filtered/filtered_deseq2.tsv"
        output:
            rds="clusterProfiler/genelist.RDS",
            translation_table="clusterProfiler/translation_table.tsv"
        message: "Translate genes identifiers to ENTREZ identifiers"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35
        params:
            gene_id=config.get("gene_id", "Gene_ID"),
            key_type=config.get("key_type", "ENSEMBL"),
            to_type=config.get("to_type", ["ENTREZID", "SYMBOL"])
        log:
            "logs/clusterProfiler/bitr_GRCh38.log"
        wrapper:
            "0.72.0-493-g8b815973b/bio/clusterProfiler/bitr_GRCh38"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/clusterProfiler/bitr_GRCh38`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

Parameters:

* gene_col    (str): Column containing gene names or identifiers
* stat_col    (str): Column containing the stat change
* cluster_col (str): Column containing the cluster to differentiate
* bdd         (str): The database to use for the analysis (default GO:BP)

Available databases:

* Gene Onthology: GO:BP, GO:CC, GO:MS
* Network of Cancer Genes: NCG
* Molecular Signature Database: MSigDB
* Disease Onthology: DO
* Disease Genes Network: DGN




Authors
-------


* Thibault

