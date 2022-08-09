.. _`Variant_Calling_Somatic`:

VARIANT_CALLING_SOMATIC
=======================

Perform Variant calling on Somatic

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_somatic/run.sh hg19

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_somatic/run.sh hg38

  # Copy/paste the following line for **MM10**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_somatic/run.sh mm10


Input/Output
------------


**Input:**

 
  
* Pairs of Tumor Fastq-formatted sequence files
  
 
  
* Pairs of corresponding Non-Tumor Fastq-formatted sequence files
  
 
  
* Cosmic database formatted as gzipped vcf and its tbi index (already provided for IGR Flamingo users)
  
 
  
* dbSNP database formatted as gzipped vcf and its tbi index (already provided for IGR Flamingo users)
  
 
  
* MSigDB database formatted as GMT (already provided for IGR Flamingo users)
  
 
  
* GWASCatalog database formatted as TSV (already provided for IGR Flamingo users)
  
 
  
* Kaviar database formatted as gzipped vcf and its tbi index (already provided for IGR Flamingo users)
  
 
  
* SnpEff database downloaded with SnpEff itself (already provided for IGR Flamingo users)
  
 
  
* CaptureKit genomic intervals formatted as BED (already provided for IGR Flamingo users)
  
 
  
* Known variants from dbSNP, with only AF within the INFO field for GATK (already provided for IGR Flamingo users)
  
 
  
* dbNSFP database formatted as TSV (already provided for IGR Flamingo users)
  
 
  
* FastQ Screen databases (already provided for IGR Flamingo users)
  
 


**Output:**

 
  
* Annotated VCF files
  
 
  
* MultiQC report
  
 





Used meta-wrappers
------------------

The following individual meta-wrappers are used in this pipeline:


* :ref:`meta/bio/bwa_fixmate`

* :ref:`meta/bio/gatk_bqsr`

* :ref:`meta/bio/varscan2_calling`

* :ref:`bigr_pipelines/snpeff_snpsift`


Please refer to each meta-wrapper in above list for additional configuration parameters and information about the executed code.




Used wrappers
-------------

The following individual wrappers are used in this pipeline:


* :ref:`bio/bigr/copy`

* :ref:`bio/fastp`

* :ref:`bio/compress/pbgzip`

* :ref:`bio/tabix`

* :ref:`bio/multiqc`

* :ref:`bio/picard/collectalignmentsummarymetrics`

* :ref:`bio/fastq_screen`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.




Notes
-----

Prerequisites for piREST:

* A TSV formatted design file, *named 'design.tsv'* with the following columns:

.. list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - Sample_id
    - Upstream_fastq
    - Downstream_fastq
    - datasetIdBED
    - genomeVersion
  * - Name of the Sample1
    - Path to upstream fastq file
    - Path to downstream fastq file
    - dataset Id corresponding to the axpected bed file
    - the genome version used for this sample
  * - Name of the Sample2
    - Path to upstream fastq file
    - Path to downstream fastq file
    - dataset Id corresponding to the axpected bed file
    - the genome version used for this sample
  * - ...
    - ...
    - ...





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    include: "rules/000.common.smk"


    rule variants_only:
        input:
            multiqc="data_output/MultiQC/Somatic_Variant_Calling.html",
            vcfs=expand("data_output/VCF/{sample}.vcf.gz", sample=sample_list),
            vcf_tbis=expand("data_output/VCF/{sample}.vcf.gz.tbi", sample=sample_list),
            vcf_xlsx=expand("data_output/XLSX/{sample}.xlsx", sample=sample_list),
            vcf_tsv=expand("data_output/TSV/{sample}.tsv", sample=sample_list),
            tmb="data_output/TMB.tsv",
            msi="data_output/MSI.tsv",
            cnv=expand("data_output/CNV/{sample}.tsv", sample=sample_list),


    rule cnv_only:
        input:
            cnv=expand("data_output/CNV/{sample}.tsv", sample=sample_list),
            multiqc="data_output/MultiQC/MappingQC.html",


    rule tmb_only:
        input:
            tmb="data_output/TMB.tsv",
            multiqc="data_output/MultiQC/Somatic_Variant_Calling.html",


    rule msi_only:
        input:
            msi="data_output/MSI.tsv",
            multiqc="data_output/MultiQC/MappingQC.html",


    rule mapping_only:
        input:
            multiqc="data_output/MultiQC/MappingQC.html",


    rule fusions_only:
        input:
            multiqc="data_output/MultiQC/FusionsQC.html",


    #####################
    ### Final archive ###
    #####################


    include: "rules/024.archive.smk"


    ###############
    ### Fusions ###
    ###############


    include: "rules/025.star_chimera.smk"
    include: "rules/026.star_fusions.smk"


    ###############################
    ### TSV and Xlsx formatting ###
    ###############################


    include: "rules/018.vcf2tsv.smk"


    #################
    ### Gather QC ###
    #################


    include: "rules/016.mapping_qc.smk"
    include: "rules/017.multiqc.smk"


    ######################
    ### MSI sensor pro ###
    ######################


    include: "rules/022.msisensor.smk"


    ###########
    ### TMB ###
    ###########


    include: "rules/021.tmb.smk"


    ##################
    ### CNV Facets ###
    ##################


    include: "rules/020.cnv_facets.smk"
    include: "rules/023.annot_sv.smk"


    ###########################
    ### VCF FILE INDEXATION ###
    ###########################


    include: "rules/019.tabix.smk"


    ######################
    ### VCF annotation ###
    ######################


    include: "rules/010.snpeff.smk"
    include: "rules/012.snpsift.smk"
    include: "rules/015.occurence.smk"
    include: "rules/011.spliceai.smk"
    include: "rules/013.vcftools.smk"
    include: "rules/014.bigr.annot.smk"


    ############################################################################
    ### Correcting Mutect2 :                                                 ###
    ### AS_FilterStatus: Number=1 and not Number=A which violates VCF format ###
    ############################################################################

    ###############################
    ### Variant calling Mutect2 ###
    ###############################


    include: "rules/009.gatk.mutect2.smk"


    ##############################
    ### GATK BAM RECALIBRATION ###
    ##############################


    include: "rules/007.gatk.recal.smk"


    #####################
    ### Deduplicating ###
    #####################


    include: "rules/006.sambamba.markdup.smk"


    # Filter a bam over the capturekit bed file


    include: "rules/005.samtools.filter.smk"


    ###################
    ### BWA MAPPING ###
    ###################


    include: "rules/003.bwa.fixmate.smk"


    ############################
    ### FASTP FASTQ CLEANING ###
    ############################


    include: "rules/002.fastp.trimming.smk"


    #################################################
    ### Gather files from iRODS or mounting point ###
    #################################################


    include: "rules/001.bigr_copy.smk"


    ###########################
    ### Datasets indexation ###
    ###########################


    include: "rules/004.sambamba.index.smk"




Authors
-------


* Thibault Dayris

* M boyba Diop

* Marc Deloger
