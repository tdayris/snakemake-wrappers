.. _`clusterp profiler`:

CLUSTERP PROFILER
=================

Run Gene set analyses with ClusterProfiler


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    default_clusterprofiler_config = {
        "organism": "Hs",
        "gsea_extra": {
            "go": "pvalueCutoff = 1, pAdjustMethod = 'BH', by = 'DOSE'"
        },
        "barplot_extra": "",
        "dotplot_extra": "",
        "heatplot_extra": "",
    }

    """
    The following rules will plot graphs on enriched/gsea objects
    """
    rule barplot:
        input:
            rds = "{method}/{db}/{comparison}/{method}.{db}.{comparison}.RDS"
        output:
            png = "results/{comparison}/barplot.{method}.{db}.png"
        message: "Barplot for {wildcards.comparison}, {wildcards.method}:{wildcards.db}"
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 15,
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
            tmpdir="tmp"
        params:
            barplot_extra=config.get("barplot_extra", "")
        log:
            "logs/barplot/{method}.{db}.{comparison}.log"
        wrapper:
            "bio/clusterProfiler/barplot"


    rule dotplot:
        input:
            rds = "{method}/{db}/{comparison}/{method}.{db}.{comparison}.RDS"
        output:
            png = "results/{comparison}/dotplot.{method}.{db}.png"
        message: "Dotplot for {wildcards.comparison}, {wildcards.method}:{wildcards.db}"
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 15,
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
            tmpdir="tmp"
        params:
            dotplot_extra=config.get("dotplot_extra", "")
        log:
            "logs/dotplot/{method}.{db}.{comparison}.log"
        wrapper:
            "bio/clusterProfiler/dotplot"

    rule heatplot:
        input:
            rds = "{method}/{db}/{comparison}/{method}.{db}.{comparison}.RDS",
            gene_list = "gene_lists/{comparison}.RDS"
        output:
            png = "results/{comparison}/heatplot.{method}.{db}.png"
        message: "Heatplot for {wildcards.comparison}, {wildcards.method}:{wildcards.db}"
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 15,
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
            tmpdir="tmp"
        params:
            heatpot_extra=config.get("heatplot_extra", "")
        log:
            "logs/heatplot/{method}.{db}.{comparison}.log"
        wrapper:
            "bio/clusterProfiler/heatplot"


    rule enrich_upsetplot:
        input:
            rds="{method}/{db}/{comparison}/{method}.{db}.{comparison}.RDS"
        output:
            png="results/{comparison}/upsetplot.{method}.{db}.png"
        message: "Enrichplot for {wildcards.comparison}, {wildcards.method}:{wildcards.db}"
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 15,
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
            tmpdir="tmp"
        params:
            upsetplot_extra=config.get("upsetplot_extra", "n = 5")
        log:
            "logs/upsetplot/{method}.{db}.{comparison}.log"
        wrapper:
            "bio/clusterProfiler/upsetplot"


    rule pathview:
        input:
            genelist="gene_lists/{comparison}.RDS"
        output:
            png = "results/{comparison}/pathview.{pid}.png"
        message: "Pathview for {wildcards.comparison}"
        threads: 1
        resources:
            time_min=lambda wildcards, attempt: attempt * 15,
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
            tmpdir="tmp"
        params:
            pathway_id=lambda wildcards: str(wildcards.pid)
        log:
            "logs/pathview/{comparison}.{pid}.log"
        wrapper:
            "bio/clusterProfiler/pathview"


    """
    The following rules will performs analysis with Gene Onthology
    """
    rule gseGO:
        input:
            rds = "gene_lists/{comparison}.RDS"
        output:
            rds = temp("gsea/{db}/{comparison}/gsea.{db}.{comparison}.RDS"),
            tsv = "results/{comparison}/gsea.{db}.{comparison}.tsv"
        message: "Running GSEA GO {params.ontology} on {wildcards.comparison}"
        threads: 1
        resources:
            time_min=lambda wildcards, input, attempt: attempt * 20,
            mem_mb=lambda wildcards, input, attempt: attempt * 1024 * 3,
            tmpdir="tmp"
        params:
            gseGO_extra = config.get(
                "gsea_extra", {"go": "pvalueCutoff = 1"}
            ).get("go"),
            organism = config.get("organism", "Hs"),
            ontology = lambda wildcards: (
                "MF" if wildcards.db == "gomf" else (
                    "CC" if wildcards.db == "gocc" else "BP"
                )
            )
        log:
            "logs/gsego/{comparison}.{db}.log"
        wrapper:
            "bio/clusterProfiler/gseGO"


    rule enrichGO:
        input:
            rds = "gene_lists/{comparison}.RDS"
        output:
            rds = temp("enrich/{db}/{comparison}/enrich.{db}.{comparison}.RDS"),
            tsv = "results/enrich.{db}.{comparison}.tsv"
        message: "Running enrichGO {params.ontology} on {wildcards.comparison}"
        threads: 1
        resources:
            time_min=lambda wildcards, input, attempt: attempt * 20,
            mem_mb=lambda wildcards, input, attempt: attempt * 1024 * 3,
            tmpdir="tmp"
        params:
            enrichGO_extra=config.get(
                "enrich_extra", {"go": ""},
            ).get("go"),
            organism = config.get("organism", "Hs"),
            ontology = lambda wildcards: (
                "MF" if wildcards.db == "gomf" else (
                    "CC" if wildcards.db == "gocc" else "BP"
                )
            )
        log:
            "logs/enrichgo/{comparison}.{db}.log"
        wrapper:
            "bio/clusterProfiler/enrichGO"


    """
    Disease Ontology analysis
    """
    rule enrichDO:
        input:
            rds = "gene_lists/{comparison}.RDS"
        output:
            rds = temp("enrich/do/{comparison}/enrich.do.{comparison}.RDS"),
            tsv = "results/enrich.do.{comparison}.tsv"
        message: "Running enrichDO on {wildcards.comparison}"
        threads: 1
        resources:
            time_min=lambda wildcards, input, attempt: attempt * 20,
            mem_mb=lambda wildcards, input, attempt: attempt * 1024 * 3,
            tmpdir="tmp"
        params:
            enrichDO_extra=config.get("enrich_extra",{"dose": ""}).get("dose", ""),
            organism = config.get("organism", "Hs"),
        log:
            "logs/enrichdo/{comparison}.log"
        wrapper:
            "bio/clusterProfiler/enrichDO"


    """
    DisGeNET analysis
    """
    rule enrichDGN:
        input:
            rds = "gene_lists/{comparison}.RDS"
        output:
            rds = temp("enrich/dgn/{comparison}/enrich.dgn.{comparison}.RDS"),
            tsv = "results/enrich.dgn.{comparison}.tsv"
        message: "Running enrichDGN on {wildcards.comparison}"
        threads: 1
        resources:
            time_min=lambda wildcards, input, attempt: attempt * 20,
            mem_mb=lambda wildcards, input, attempt: attempt * 1024 * 3,
            tmpdir="tmp"
        params:
            enrichDGN_extra=config.get("enrich_extra",{"dgn": ""}).get("dgn", ""),
            organism = config.get("organism", "Hs"),
        log:
            "logs/enrichdgn/{comparison}.log"
        wrapper:
            "bio/clusterProfiler/enrichDGN"


    """
    Network of Cancer Genes analysis
    """
    rule enrichNCG:
        input:
            rds = "gene_lists/{comparison}.RDS"
        output:
            rds = temp("enrich/ncg/{comparison}/enrich.ncg.{comparison}.RDS"),
            tsv = "results/enrich.ncg.{comparison}.tsv"
        message: "Running enrichNCG on {wildcards.comparison}"
        threads: 1
        resources:
            time_min=lambda wildcards, input, attempt: attempt * 20,
            mem_mb=lambda wildcards, input, attempt: attempt * 1024 * 3,
            tmpdir="tmp"
        params:
            enrichDGN_extra=config.get("enrich_extra",{"ncg": ""}).get("ncg", ""),
            organism = config.get("organism", "Hs"),
        log:
            "logs/enrichncg/{comparison}.log"
        wrapper:
            "bio/clusterProfiler/enrichNCG"


    """
    Molecular Signature database analysis
    """
    rule gsea_MSigDB:
        input:
            rds = "gene_lists/{comparison}.RDS"
        output:
            rds = temp("gsea/msigdb_{category}/{comparison}/gsea.msigdb_{category}.{comparison}.RDS"),
            tsv = "results/gsea.msigdb_{category}.{comparison}.tsv"
        message: "Running MSigDB {wildcards.category} on {wildcards.comparison}"
        threads: 1
        resources:
            time_min=lambda wildcards, input, attempt: attempt * 20,
            mem_mb=lambda wildcards, input, attempt: attempt * 1024 * 3,
            tmpdir="tmp"
        params:
            msigdb_extra = lambda wildcards: f"category={wildcards.category}",
            organism = config.get("organism", "Hs"),
            msigdb_gsea_extra = config.get("gsea_extra",{"msigdb": ""}).get("msigdb", ""),
        log:
            "logs/gseamsigdb/{comparison}_{category}.log"
        wrapper:
            "bio/clusterProfiler/enrichNCG"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/clusterProfiler/enrichDGN`

* :ref:`bio/clusterProfiler/enrichNCG`

* :ref:`bio/clusterProfiler/enrichDO`

* :ref:`bio/clusterProfiler/enrichGO`

* :ref:`bio/clusterProfiler/gseGO`

* :ref:`bio/clusterProfiler/pathview`

* :ref:`bio/clusterProfiler/upsetplot`

* :ref:`bio/clusterProfiler/heatplot`

* :ref:`bio/clusterProfiler/dotplot`

* :ref:`bio/clusterProfiler/barplot`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.







Authors
-------


* Thibault Dayris

