.. _`sarek enhanced`:

SAREK ENHANCED
==============

Run Sarek and add supplementary annotations, better default parameters, automatic fastq detection, automatic fastq retrivement from iRODS and GLeaves complience.

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 


Input/Output
------------


**Input:**

 
  
* Fastq files (optional, if not provided, then a file called `design.tsv` must contain a path to the fastq files)
  
 
  
* TSV formatted design file (optional, if not provided, the fastq files must be available in the current working directory tree)
  
 


**Output:**

 
  
* Sarek results with enhancements
  
 









Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    localrules: sarek

    # TODO: Create input file
    # TODO: Create flamingo configuration


    rule sarek:
        input:
            nf=config["nextflow"]["nf"],
            input="input.tsv",
            igenomes_base=config["sarek"]["db"],
            flamingo_config=config["nextflow"]["flamingo"]
        output:
            directory("work"),
            directory("results"),
            temp(".nextflow_log")
        message:
            "Running NF-core Sarek 2.7.1"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            config["sarek"].get("extra", "")
        handover: True
        envmodules:
            'nextflow/21.04.3',
            'singularity/3.6.3',
            'java/12.0.2'
        log:
            "logs/nextflow/sarek.log"
        shell:
            "declare -x NXF_OPTS='-Xms1g -Xmx4g' && "
            "export NXF_OPTS && "
            "nextflow run {input.nf} "
            "--input {input.input} "
            "--igenomes_base {input.igenomes_base} "
            "{extra} "
            "--max_cpus {threads} "
            "--max_memory {resources.mem_mb} "
            "--max_time {resources.time_min}"
            "> {log} 2>&1"




Authors
-------


* Thibault Dayris
