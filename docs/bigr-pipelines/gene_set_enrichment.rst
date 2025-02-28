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

    worflow_source_dir = Path(snakemake.workflow.srcdir("."))
    common = str(worflow_source_dir / ".." / "common" / "python")
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

    dbs = list(config.get("gmt", {}).keys()) + [
        #"NetworkCancerGenes",
        #"DiseaseOnt",
        #"DisGenNet"
    ]

    ppis = list(config.get("ppi", {}).keys())
    all_dbs = dbs + ppis

    methods = ["enrich"] # ["gsea", "enrich"]
    plots = ["barplot", "dotplot", "upsetplot"] # "heatplot", "treeplot"]
    keytypes = ["ENSEMBL", "ENTREZID", "SYMBOL", "ENSEMBLPROT"]

    db_key = {}
    for db in dbs:
        if db in config["gmt"].keys():
            if config["gmt"][db].endswith(".entrez.gmt"):
                    db_key[db] = "ENTREZID"
                elif config["gmt"][db].endswith(".ENSG.gmt"):
                    db_key[db] = "ENSEMBL"
                elif config["gmt"][db].endswith(".symbols.gmt"):
                    db_key[db] = "SYMBOL"
            else:
                db_key[db] = "ENTREZID"


    def plot_list():
        results = []
        key = None

        for db in dbs:
            results.append(
                f"results/all_comparisons/{db}.{db_key[db]}.png"
            )
            for met in methods:
                for plot in plots:
                    for comp in comparisons:
                        results.append(
                            f"results/{comp}/{db}.{db_key[db]}/{plot}.{met}.png"
                        )

        for ppi in ppis:
            results.append(
                f"results/all_comparisons/{db}.ENSEMBLPROT.png"
            )
            for met in methods:
                for plot in plots:
                    for comp in comparisons:
                        results.append(f"results/{comp}/{ppi}.ENSEMBLPROT/{plot}.{met}.png")

        return results

    def tsv_list():
        results = []
        for comparison in comparisons:
            for db in dbs:
                results.append(
                    f"results/{comparison}/{db}.{db_key[db]}/enrich.{comparison}.{db_key[db]}.tsv"
                )


    ruleorder: enrichDO > enricherGMT
    ruleorder: enrichDGN > enricherGMT
    ruleorder: enrichNCG > enricherGMT


    wildcard_constraints:
        comparison=r"|".join(comparisons),
        method=r"|".join(methods),
        db=r"|".join(all_dbs),
        keytype=r"|".join(keytypes),
        ppi=r"|".join(ppis)


    rule target:
        input:
            plot_list()



    rule expand_rank_list:
        input:
            tsv=config.get("ranks.list.tsv", "ranks.list.tsv")
        output:
            tsv=expand(
                "gene_lists/ENTREZID/{comparison}.tsv",
                comparison=comparisons
            ),
            entrez_rds=temp(expand(
                "gene_lists/ENTREZID/{comparison}.RDS",
                comparison=comparisons
            )),
            symbol_rds=temp(expand(
                "gene_lists/SYMBOL/{comparison}.RDS",
                comparison=comparisons
            )),
            ensembl_rds=temp(expand(
                "gene_lists/ENSEMBL/{comparison}.RDS",
                comparison=comparisons
            )),
            protein_rds = temp(expand(
                "gene_lists/ENSEMBLPROT/{comparison}.RDS",
                comparison=comparisons
            )),
            universe = temp(expand(
                "gene_lists/universe/{comparison}.RDS",
                comparison=comparisons
            ))
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


    ###########################################
    ### Plot comparisons agains each others ###
    ###########################################

    rule paired_dotplots:
        input:
            tsv = expand(
                "results/{comparison}/{database}.{keytype}/enrich.{comparison}.{keytype}.tsv",
                comparison=comparisons,
                allow_missing=True
            )
        output:
            png = "results/all_comparisons/{database}.{keytype}.png"
        message: "Aggregating multiple comparison together"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/paired_dotplot/{database}.{keytype}.log"
        conda:
            "envs/seaborn.yaml"
        params:
            comparisons=comparisons
        script:
            "scripts/relplot.py"


    #################
    ### Big table ###
    #################

    rule concat_tables:
        input:
            tsv = tsv_list()
        output:
            temp("tmp/complete_table_wh.tsv")
        message:
            "Joining all TSV results"
        threads: 1
        resources:
            mem_mb=256,
            time_min=lambda wildcards, attempt * 5,
            tmpdir="tmp"
        log:
            "logs/awk/concat.log"
        params:
            begin=['FS=OFS="\t"'],
            body=["{if (NR==1) {print $0 FS FileName} else {print $0 FS FILENAME}}"]
        wrapper:
            "bio/awk"


    rule clean_headers:
        input:
            "tmp/complete_table_wh.tsv"
        output:
            "results/complete_table.tsv"
        message:
            "Cleaning remaining headers in joint table"
        threads: 1
        resources:
            mem_mb=256,
            time_min=lambda wildcards, attempt * 5,
            tmpdir="tmp"
        log:
            "logs/sed/concat.log"
        params:
            regex="1p;s/^ID\tDescri//g"
        shell:
            "sed {params.regex} {input} > {output} 2> {log}"






Authors
-------


* Thibault Dayris
