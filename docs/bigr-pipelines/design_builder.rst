.. _`design_builder`:

DESIGN_BUILDER
==============

Build design file from sample list

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 


Input/Output
------------


**Input:**

 
  
* A list of samples, one per line, stored in a file called "samples.txt".
  
 


**Output:**

 
  
* A TSV-formatted design file with Sample_id, Upstream_file, Downstream_file
  
 






Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`iRODS/search_samples`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.




Notes
-----

The samples must be in iRODS.

Under development. Do no use.





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    rule name_cols:
        input:
            "irods/filtered.tsv"
        output:
            "design.tsv"
        group:
            "format"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 512,
            time_min=lambda wildcards, attempt: attempt * 5
        log:
            "logs/irods/cols.log"
        params:
            ""
        shell:
            'cat <(echo -e "Upstream_file\tSample_id\tDownstream_file") {input} > {output} 2> {log}'


    rule filter_results:
        input:
            table="irods/table.tsv",
            samples="samples.txt"
        output:
            temp("irods/filtered.tsv")
        group:
            "format"
        threads: 7
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 512,
            time_min=lambda wildcards, attempt: attempt * 5
        log:
            "logs/irods/filter.log"
        params:
            re='mRNA-seq',
            subset='-f3,16'
        shell:
            "grep --file={input.samples} {input.table} | grep {params.re} | cut {params.subset} | sort | uniq | paste - - | cut -f1-3 > {output} 2> {log}"


    rule search_samples:
        input:
            samples="samples.txt"
        output:
            "irods/table.tsv"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 25
        log:
            "logs/irods/search_samples.log"
        wrapper:
            "bio/iRODS/search_samples"




Authors
-------


* Thibault Dayris
