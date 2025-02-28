.. _`nfcore-chipseq`:

NFCORE-CHIPSEQ
==============

Perform various analyses on Chipseq

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your directory containing VCF files

  cd /path/to/project

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/nfcore_chipseq/run.sh


Input/Output
------------


**Input:**

 
  
* Fastq files
  
 


**Output:**

 
  
* Peak called
  
 









Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    from snakemake.utils import min_version
    from pathlib import Path
    from yaml import dump
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
        filename="snakemake.salmon_quant.log",
        filemode="w",
        level=logging.DEBUG
    )

    default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
    configfile: get_config(default_config)
    design = get_design(os.getcwd(), search_fastq_pairs)

    #try:
    fastq_links = link_fq(
        design.Sample_id,
        design.Upstream_file,
        design.Downstream_file
    )

    default_config = """//Profile config names for nf-core/configs
    params {
      config_profile_description = 'The Flamingo cluster profile'
      config_profile_contact = 'Thibault Dayris (@BiGR)'
      config_profile_url = 'https://gitlab.com/bioinfo_gustaveroussy/bigr'
    }

    singularity {
      enabled = true
      autoMounts = false
      runOptions = '-B /mnt/beegfs:/mnt/beegfs'
    }

    process {
      executor = 'slurm'
    }

    params {
      igenomes_ignore = true
      igenomesIgnore = true //deprecated
      max_memory = 750.GB
      max_cpus = 200
      max_time = 24.h
    }"""


    onstart:
        shell("export NXF_SINGULARITY_CACHEDIR=/mnt/beegfs/software/nf-core-chipseq/1.2.2/")





    include: "rules/000.get_fastq.smk"
    include: "rules/001.prepare_nfcore.smk"
    include: "rules/002.nfcore_chipseq.smk"

    rule target:
        input:
            "results/mulitqc/broadPeak/multiqc_report.html"




Authors
-------


* Thibault Dayris
