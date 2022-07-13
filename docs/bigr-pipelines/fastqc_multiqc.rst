.. _`FastQC_MultiQC`:

FASTQC_MULTIQC
==============

Simply run FastQC and MultiQC on available Fastq files

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to the directory containing your fastq files (they may be in sub-directories)

  cd /path/to/my/fastq.gz

  # Copy paste the following line

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/fastqc_multiqc/run.sh


Input/Output
------------


**Input:**

 
  
* Fastq files
  
 
  
* FastQ-Screen databases (already provided for IGR Flamingo users)
  
 


**Output:**

 
  
* HTML report for each single fastq files
  
 
  
* Complete HTML report for all fastq files
  
 






Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`bio/bigr/copy`

* :ref:`bio/fastqc`

* :ref:`bio/multiqc`

* :ref:`bio/fastq_screen`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.




Notes
-----

This simple pipeline is here to assess the quality after automatic sequencing process. Keep in mind that FastQC is built for DNASeq, and may raise non-critical warnings.

This pipeline requires either (1) a two columned design file called 'design.tsv', or (2) available fastq files in the current directory.

.. list-table:: Desgin file format
    :widths: 33 33
    :header-rows: 1

    * - Sample_id
      - Upstream_fastq
    * - Name of the Sample1
      - Path to upstream fastq file
    * - Name of the Sample2
      - Path to upstream fastq file
    * - ...
      - ...





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    import sys
    from pathlib import Path

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../common/python")
    sys.path.append(common)

    from dataframes import *
    from file_manager import *
    from files_linker import *
    from graphics import *
    from write_yaml import *
    from messages import message

    from snakemake.utils import min_version
    min_version("6.0")

    logging.basicConfig(
        filename="snakemake.fastqc_multiqc.log",
        filemode="w",
        level=logging.DEBUG
    )

    default_config = read_yaml(worflow_source_dir / "config.yaml")
    config_path = get_config(default_config)
    design = get_design(os.getcwd(), search_fastq_pairs)


    fastq_links = link_fq(
        design.Sample_id,
        design.Upstream_file,
        design.Downstream_file
    )

    configfile: config_path
    container: "docker://continuumio/miniconda3:4.4.10"


    ##################################
    ### Gather all quality reports ###
    ##################################

    rule multiqc:
        input:
            fqc_zip=expand(
                "fastqc/{sample}_{stream}_fastqc.zip",
                sample=design["Sample_id"],
                stream=["1", "2"]
            ),
            fqc_html=expand(
                "fastqc/{sample}.{stream}.html",
                sample=design["Sample_id"],
                stream=["1", "2"]
            ),
            txt=expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.txt",
                sample=design["Sample_id"],
                stream=["1", "2"]
            ),
            png=expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.png",
                sample=design["Sample_id"],
                stream=["1", "2"]
            )
        output:
            "multiqc/multiqc.html",
            directory("multiqc/multiqc_data")
        message:
            "Gathering all quality reports in {output}"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: attempt * 2048,
            time_min=lambda wildcard, attempt: attempt * 50,
            tmpdir="tmp"
        params:
            "--flat"
        log:
            "logs/multiqc.log"
        wrapper:
            "bio/multiqc"


    #######################################################
    ### Adding specific actions for BiGR demultiplexing ###
    #######################################################

    use rule multiqc as irods_complient with:
        input:
            fqc_zip=expand(
                "fastqc/{sample}_{stream}_fastqc.zip",
                sample=design["Sample_id"],
                stream=["1", "2"]
            ),
            fqc_html=expand(
                "fastqc/{sample}.{stream}.html",
                sample=design["Sample_id"],
                stream=["1", "2"]
            ),
            txt=expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.txt",
                sample=design["Sample_id"],
                stream=["1", "2"]
            ),
            png=expand(
                "fastq_screen/{sample}.{stream}.fastq_screen.png",
                sample=design["Sample_id"],
                stream=["1", "2"]
            ),
            bcl_json="Stats.json"
        output:
            "output/multiqc.html",
            directory("output/multiqc_data")
        group:
            "stats_inclusion"


    rule unzip_stats:
        output:
            temp("Stats.json")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/unzipping.log"
        group:
            "stats_inclusion"
        shell:
            'unzip -n -d "${{PWD}}" '
            'input/*/archive/*/unaligned/Stats/Stats.json.zip '
            '> {log} 2>&1'


    #########################################
    ### Assess quality of each fastq file ###
    #########################################

    rule fastqc:
        input:
            "reads/{sample}.{stream}.fq.gz"
        output:
            html=temp("fastqc/{sample}.{stream}.html"),
            zip=temp("fastqc/{sample}_{stream}_fastqc.zip")
        message:
            "Assessing quality of {wildcards.sample}, ({wildcards.stream})"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 4096),
            time_min=lambda wildcard, attempt: attempt * 50,
            tmpdir="tmp"
        params:
            ""
        log:
            "logs/fastqc/{sample}.{stream}.log"
        wrapper:
            "bio/fastqc"


    rule fastq_screen:
        input:
            "reads/{sample}.{stream}.fq.gz"
        output:
            txt=temp("fastq_screen/{sample}.{stream}.fastq_screen.txt"),
            png=temp("fastq_screen/{sample}.{stream}.fastq_screen.png")
        message:
            "Assessing quality of {wildcards.sample}, stream {wildcards.stream}"
        threads: config.get("threads", 20)
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 8192),
            time_min=lambda wildcard, attempt: attempt * 50,
            tmpdir="tmp"
        params:
            fastq_screen_config=config["fastq_screen"],
            subset=100000,
            aligner='bowtie2'
        log:
            "logs/fastq_screen/{sample}.{stream}.log"
        wrapper:
            "bio/fastq_screen"


    #################################################
    ### Gather files from iRODS or mounting point ###
    #################################################

    rule bigr_copy:
        output:
            "reads/{sample}.fq.gz"
        message:
            "Gathering {wildcards.sample} fastq file"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
            time_min=lambda wildcard, attempt: attempt * 45,
            tmpdir="tmp"
        params:
            input=lambda wildcards, output: fastq_links[output[0]]
        log:
            "logs/bigr_copy/{sample}.log"
        wrapper:
            "bio/BiGR/copy"




Authors
-------


* Thibault Dayris

* Gérôme Jules-Clément

* Marie Martelat
