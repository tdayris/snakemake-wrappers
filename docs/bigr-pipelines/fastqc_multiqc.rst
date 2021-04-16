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
    :widths: 33 33 33
    :header-rows: 1

    * - Sample_id
      - Upstream_fastq
    * - Name of the Sample1
      - Path to upstream fastq file
    * - Name of the Sample2
      - Path to upstream fastq file
    * - ...
      - ...




Authors
-------


* Thibault Dayris

* Gérôme Jules-Clément

* Marie Martelat
