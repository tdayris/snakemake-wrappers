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


How does it work ?

1. Gathering Fastq files

This pipeline copies files from iRODS, or symlinks them. In your design file,
the `Sample_id` is used to rename your fastq file and make them all easily
recognisable. In case of iRODS copy, the checksum is automatically computed
and verified. In cas of a copy from cold to hot storage, then a checksum is
automatically computed and verified. Elsewise, a simple sylmink is done and
no control needs to be performed.

The column Upstream_file/Downstream_file identifies reads' streams.
If the sequencing was not oriented, then order does not matter.
Otherwise, make sure R1 reads are under Upstream_file, and R2 reads under
Downstream_file.

You may need to concatenate several fastq files into one single fastq file
for a given sample: in case of lane splitting, run splitting, and/or
resequencing. This may be done automatically! Under the corresponding column,
separate the multiple files by a comma (`,`).

2. Cleaning fastq files

These (concatenated?) fastq files are trimmed with fastp. By default, the
cleaning is a bit more relaxed than default fastp parameters. We're woring on
RNA-Seq with Salmon in this pipeline. No need to be very strict.

The running window has been increased to 6 nucleotides, the minimum mean read
quality was raised to 10. The unqualified percent limit was raised to 50% and
the maximum of N bases was raised to 7. The minimum length was lowered to 15.
Because we are trimming RNA-Seq, we analyse overrepresented sequences.

Nothing dramatic! Bad reads will be filtered later if needed.

By the way, FastQ Screen is used over the raw reads, since possible sequencing
artifacts are interesting in this step. Many genomes are tested, see section
5 to look at the results.

3. Gathering quality reports

MultiQC aggregates all quality reports so you can compare all your samples
easily.





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    import os
    import sys
    from pathlib import Path
    from snakemake.utils import min_version
    min_version("7.0")


    worflow_source_dir = Path(snakemake.workflow.srcdir("."))
    common = str(worflow_source_dir / ".." / "common" / "python")
    sys.path.append(common)


    from dataframes import *
    from file_manager import *
    from files_linker import *
    from graphics import *
    from write_yaml import *
    from reservation import *
    from messages import message

    from snakemake.utils import min_version
    from snakemake.shell import shell
    min_version("6.0")

    logging.basicConfig(
        filename="snakemake.fastqc_multiqc.log",
        filemode="w",
        level=logging.DEBUG
    )

    default_config = read_yaml(worflow_source_dir / "config.yaml")
    config_path = get_config(default_config)
    design = get_design(os.getcwd(), search_fastq_files)
    try:
        design.columns = ["Sample_id", "Upstream_file"]
    except ValueError:
        pass

    fastq_links = link_fq(
        design.Sample_id,
        design.Upstream_file,
    )

    configfile: config_path
    container: "docker://continuumio/miniconda3:4.4.10"


    ##################
    ### Flag rules ###
    ##################

    onsuccess:
        shell("touch DONE && rm --force --verbose ON_GOING ERROR")

    onerror:
        shell("touch ERROR && rm --force --verbose ON_GOING DONE")

    onstart:
        shell("touch ON_GOING && rm --force --verbose ERROR DONE")


    ##################################
    ### Gather all quality reports ###
    ##################################


    include: "rules/003.multiqc.smk"


    #########################################
    ### Assess quality of each fastq file ###
    #########################################


    include: "rules/001.fastqc.smk"
    include: "rules/002.fastq_screen.smk"


    #################################################
    ### Gather files from iRODS or mounting point ###
    #################################################


    include: "rules/000.copy.smk"


    prefix = "multiqc"


    rule target:
        input:
            f"{prefix}/multiqc.html"




Authors
-------


* Thibault Dayris

* Gérôme Jules-Clément

* Marie Martelat
