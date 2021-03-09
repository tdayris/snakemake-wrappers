.. _`cluster_profiler`:

CLUSTER_PROFILER
================

Perform gene sets analysis against pathways


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    rule translate_genes_identifiers_hg38:
        input:
            table = "deseq2/filtered/filtered_deseq2.tsv"
        output:
            rds = "clusterProfiler/genelist.RDS",
            translation_table = "clusterProfiler/translation_table.tsv"
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

