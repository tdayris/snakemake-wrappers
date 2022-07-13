.. _`gene_set_enrichment (Under development)`:

GENE_SET_ENRICHMENT (UNDER DEVELOPMENT)
=======================================

Perform gene set enrichment analysis and produce wide range of graphs

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 


Input/Output
------------


**Input:**

 
  
* A TSV-formatted text table with one line per gene, one column per factor of interest and a numerical value in each of them.
  
 
  
* An optional text file containing list of pathways of interest. Pathway names must be exact.
  
 


**Output:**

 
  
* Gene Set enrichment table
  
 
  
* multiple graphs
  
 







Notes
-----

This pipeline can be used in "two times". First run the pipeline from the TSV input file, then look at the pathways of interest, then relaunch the pipeline to include focuced graphs and comparisons.

If only one comparison is provided in the input TSV file, then a single Gene set enrichment analysis is done. If more than one comparison is done, then each of them gets a separated gene set enrichment analysis.

The input TSV file looks like:

.. list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - Ensembl_Gene_ID
    - Comparinson1
    - Optional comparison 2
    - ...
  * - ENSX000000XXXXXX
    - FoldChange/PAdj
    - FoldChange/PAdj
    - ...
  * - ENSX000000YYYYYY
    - FoldChange/PAdj
    - FoldChange/PAdj
    - ...
  * - ...
    - ...
    - ...
    - ...





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    from snakemake.utils import min_version
    from pathlib import Path
    from yaml import dump
    from typing import List
    min_version("6.0")

    import sys

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../common/python")
    sys.path.append(common)

    from file_manager import *
    from files_linker import *
    from write_yaml import *
    from messages import message

    logging.basicConfig(
        filename="snakemake.gene_set_enrichment.log",
        filemode="w",
        level=logging.DEBUG
    )

    default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
    configfile: get_config(default_config)


    def get_comparison_list(path: str) -> List[str]:
        """
        Return list of comparison levels
        Open the input file listing comparisons and ranking value
        """
        with open(path, "r") as ranks:
            header = next(ranks)
            return header.strip().split("\t")[1:]

    comparisons = get_comparison_list(
        config.get("ranks.list.tsv", "ranks.list.tsv")
    )

    dbs = ["gomf", "ncg", "gocc", "gobp", "do", "dgn"]
    methods = ["gsea", "enrich"]
    plots = ["barplot", "dotplot", "upsetplot", "heatplot"]


    ruleorder: enrichDO > enrichGO
    ruleorder: enrichDGN > enrichGO
    ruleorder: enrichNCG > enrichGO

    wildcard_constraints:
        comparison=r"|".join(comparisons),
        method=r"|".join(methods),
        db=r"|".join(dbs)


    rule target:
        input:
            rds=expand(
                "gene_lists/{comparison}.RDS",
                comparison=comparisons
            ),
            png_enrich=expand(
                "results/{comparison}/{plot}.{method}.{db}.png",
                comparison=comparisons, plot=plots, method=methods[1], db=dbs[2:]
            ),
            # png_go=expand(
            #     "results/{comparison}/{plot}.{method}.{db}.png",
            #     comparison=comparisons, plot=plots[2], method=methods[0], db=dbs
            # )



    rule expand_rank_list:
        input:
            tsv=config.get("ranks.list.tsv", "ranks.list.tsv")
        output:
            tsv=expand(
                "gene_lists/{comparison}.tsv",
                comparison=comparisons
            ),
            rds=expand(
                "gene_lists/{comparison}.RDS",
                comparison=comparisons
            )
        message: "Expanding rank lists"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
            time_min=lambda wildcards, attempt: attempt * 20,
            tmpdir="tmp"
        params:
            gene_id_type=config.get(
                "gene_id_type", "ENSEMBL"
            )
        log:
            "logs/expand.log"
        wrapper:
            "bio/clusterProfiler/hg38_genelist"


    ##############################
    ### Use of ClusterProfiler ###
    ##############################


    module clusterprofiler_meta:
        snakefile: "../../meta/bio/clusterprofiler/test/Snakefile"
        config: config


    use rule * from clusterprofiler_meta




Authors
-------


* Thibault Dayris
