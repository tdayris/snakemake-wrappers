.. _`Variant_Calling_Somatic (under development)`:

VARIANT_CALLING_SOMATIC (UNDER DEVELOPMENT)
===========================================

Perform Variant calling on Somatic

Usage
-----

In order to run the pipeline, use the following commands

.. code-block:: bash 

  # Go to your working directory

  cd /path/to/my/working/directory

  # Build a design file (see below)

  # Copy/paste the following line for **HG19**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_somatic/run.sh

  # Copy/paste the following line for **HG38**

  bash /mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/variant_calling_somatic/run.sh hg38


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

Prerequisites:

* A TSV formatted design file, *named 'design.tsv'* with the following columns:

.. list-table:: Desgin file format
  :widths: 33 33 33
  :header-rows: 1

  * - Sample_id
    - Upstream_fastq
    - Downstream_fastq
  * - Name of the Sample1
    - Path to upstream fastq file
    - Path to downstream fastq file
  * - Name of the Sample2
    - Path to upstream fastq file
    - Path to downstream fastq file
  * - ...
    - ...
    - ...





Snakefile
---------

The pipeline contains the following steps:

.. code-block:: python

    import logging
    import os
    import pandas
    import sys
    from pathlib import Path

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../common/python")
    sys.path.append(common)

    from file_manager import *
    from files_linker import *
    from write_yaml import *
    from messages import *
    from snakemake.utils import min_version
    min_version("6.0")

    logging.basicConfig(
        filename="snakemake.variant_calling_somatic.log",
        filemode="w",
        level=logging.DEBUG
    )

    container: "docker://continuumio/miniconda3:4.4.10"
    localrules: bigr_copy

    ruleorder: samtools_index_bam > samtools_index

    default_config = read_yaml(worflow_source_dir / "config.hg38.yaml")
    configfile: get_config(default_config)
    design = get_design(os.getcwd(), search_fastq_pairs)
    # design = design.head(2).tail(1)
    design.dropna(inplace=True)
    print(design)

    wildcard_constraints:
        sample = r"|".join(design["Sample_id"]),
        stream = r"1|2|R1|R2",
        status = r"normal|tumor",
        content = r"snp|indel"


    fastq_links = link_fq_somatic(
        sample_names=design.Sample_id,
        n1_paths=design.Upstream_file_normal,
        t1_paths=design.Upstream_file,
        n2_paths=design.Downstrea_file_normal,
        t2_paths=design.Downstream_file,
    )


    rule all:
        input:
            # calls=expand(
            #     "snpsift/dbnsfp/{sample}.vcf.gz{index}",
            #     sample=design["Sample_id"].tolist(),
            #     index=["", ".tbi"]
            # ),
            # html="multiqc/variant_calling_somatic.html",
            mapped=expand(
                "picard/markduplicates/{sample}_{status}.bam{ext}",
                sample=design["Sample_id"].tolist(),
                status=["normal", "tumor"],
                ext=["", ".bai"]
            ),
            # mutect2=expand(
            #     "mutect2/filter/{sample}.vcf.gz",
            #     sample=design["Sample_id"].tolist()
            # ),
            varscan2=expand(
                "varscan2/concat/{sample}.vcf.gz",
                sample=design["Sample_id"].tolist()
            ),
            qc="multiqc/variant_calling_somatic.html"
        message:
            "Finishing the WES Somatic Variant Calling"


    #################
    ### Gather QC ###
    #################

    rule multiqc:
        input:
            html=expand(
                "fastp/html/pe/{sample}_{status}.fastp.html",
                sample=design["Sample_id"],
                status=["normal", "tumor"]
            ),
            json=expand(
                "fastp/json/pe/{sample}_{status}.fastp.json",
                sample=design["Sample_id"],
                status=["normal", "tumor"]
            ),
            picard_metrics=expand(
                "picard/metrics/{sample}_{status}.picard.metrics.txt",
                sample=design["Sample_id"],
                status=["normal", "tumor"]
            ),
            fastq_screen=expand(
                "fastq_screen/{sample}.{stream}.{status}.fastq_screen.{ext}",
                sample=design["Sample_id"],
                stream=["1", "2"],
                ext=["txt", "png"],
                status=["normal", "tumor"]
            ),
            picard_summary=expand(
                "picard/alignment_summary/{sample}_{status}.summary.txt",
                sample=design["Sample_id"],
                status=["normal", "tumor"]
            )
        output:
            report(
                "multiqc/variant_calling_somatic.html",
                caption="../common/reports/multiqc.rst",
                category="Quality Controls"
            )
        message:
            "Aggregating quality reports from SnpEff"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
            time_min=lambda wildcards, attempt: attempt * 35,
            tmpdir="tmp"
        log:
            "logs/multiqc.log"
        wrapper:
            "bio/multiqc"


    rule alignment_summary:
        input:
            bam="samtools/sort/{sample}_{status}.bam",
            bam_index=get_bai("samtools/sort/{sample}_{status}.bam"),
            ref=config['ref']['fasta'],
            ref_idx=get_fai(config['ref']['fasta']),
            ref_dict=get_dict(config['ref']['fasta']),
        output:
            temp("picard/alignment_summary/{sample}_{status}.summary.txt")
        message:
            "Collecting alignment metrics on GATK recalibrated {wildcards.sample}"
            " (considering {wildcards.status})"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1020,
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir="tmp"
        log:
            "logs/picard/alignment_summary/{sample}_{status}.log"
        params:
            "VALIDATION_STRINGENCY=LENIENT "
            "METRIC_ACCUMULATION_LEVEL=null "
            "METRIC_ACCUMULATION_LEVEL=SAMPLE"
        wrapper:
            "bio/picard/collectalignmentsummarymetrics"


    rule fastq_screen:
        input:
            "reads/{status}/{sample}.{stream}.fq.gz"
        output:
            txt=temp("fastq_screen/{sample}.{stream}.{status}.fastq_screen.txt"),
            png=temp("fastq_screen/{sample}.{stream}.{status}.fastq_screen.png")
        message:
            "Assessing quality of {wildcards.sample}, {wildcards.stream}"
            " (considering {wildcards.status})"
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
            "logs/fastqc/{sample}.{stream}.{status}.log"
        wrapper:
            "bio/fastq_screen"


    ##################
    ### CNV Facets ###
    ##################


    rule cnv_facets:
        input:
            tumor_bam="picard/markduplicates/{sample}_tumor.bam",
            tumor_bai=get_bai("picard/markduplicates/{sample}_tumor.bam"),
            normal_bam="picard/markduplicates/{sample}_normal.bam",
            normal_bai=get_bai("picard/markduplicates/{sample}_normal.bam"),
            vcf=config["ref"]["dbsnp"],
            vcf_index=get_tbi(config["ref"]["dbsnp"]),
            bed=config["ref"]["capture_kit_bed"]
        output:
            vcf="facets/{sample}/{sample}.vcf.gz",
            profile="facets/{sample}/{sample}.cnv.png",
            coverage="facets/{sample}/{sample}.cov.pdf",
            spider="facets/{sample}/{sample}.spider.pdf",
            pileup="facets/{sample}/{sample}.csv.gz"
        message:
            "Searching for CNV in {wildcards.sample} with Facets"
        threads: 20
        resources:
          mem_mb=lambda wildcards, attempt: attempt * 8192,
          time_min=lambda wildcards, attempt: attempt * 30,
          tmpdir="tmp"
        params:
            extra="--snp-count-orphans"
        log:
            "logs/facets/cnv/{sample}.log"
        wrapper:
            "0.77.0-840-gff109e8d3/bio/facets/cnv"


    #################################
    ### FINAL VCF FILE INDEXATION ###
    #################################

    module compress_index_vcf_meta:
        snakefile: "../../meta/bio/compress_index_vcf/test/Snakefile"
        config: config

    use rule * from compress_index_vcf_meta


    use rule tabix_index from compress_index_vcf_meta as snp_indel_tabix_index with:
        input:
            "{tool}/{subcommand}/{sample}.{content}.vcf.gz"
        output:
            "{tool}/{subcommand}/{sample}.{content}.vcf.gz.tbi"
        message:
            "Indexing {wildcards.sample} (from somatic varscan "
            "{wildcards.content}) with tabix."
        log:
            "logs/{tool}/{subcommand}/tabix/index/{sample}.{content}.log"


    use rule pbgzip_compress from compress_index_vcf_meta as si_pbgzip with:
        input:
            "{tool}/{subcommand}/{sample}.{content}.vcf"
        output:
            "{tool}/{subcommand}/{sample}.{content}.vcf.gz"
        message:
            "Compressnig {wildcards.sample} (from somatic varscan "
            "{wildcards.content}) with pbgzip."
        log:
            "logs/{tool}/{subcommand}/pgbzip/varcsanc2/{sample}.{content}.log"


    ######################
    ### VCF annotation ###
    ######################


    module snpeff_meta:
        snakefile: "../../meta/bio/snpeff_annotate/test/Snakefile"
        config: config

    use rule snpeff from snpeff_meta with:
        input:
            calls="meta_caller/calls/{sample}.vcf.gz",
            calls_index=get_tbi("meta_caller/calls/{sample}.vcf.gz"),
            db=config["ref"]["snpeff"]


    module snpsift:
        snakefile: "../../meta/bio/snpsift/test/Snakefile"
        config: config


    use rule * from snpsift as *


    #####################################
    ### Merge variant calling results ###
    #####################################

    # module metacaller_somatic_meta:
    #     snakefile: "../../meta/bio/meta_caller_somatic/test/Snakefile"
    #     config: {"genome": config["ref"]["fasta"], "bed": config["ref"]["capture_kit_bed"]}
    #
    #
    # use rule * from metacaller_somatic_meta as *


    ############################################################################
    ### Correcting Mutect2 :                                                 ###
    ### AS_FilterStatus: Number=1 and not Number=A which violates VCF format ###
    ### AD becomes ADM: AD is reserved for Allele Depth, Mutect2 stores      ###
    ###                 multiple information under "AD" field.               ###
    ############################################################################

    rule correct_mutect2_vcf:
        input:
            "mutect2/filter_reheaded/{sample}.vcf.gz"
        output:
            temp("mutect2/corrected/{sample}.vcf")
        message:
            "Renaming reserved AD field and fixing AS_FilterStrand format error"
            " on {wildcards.sample}"
        threads: 3
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 256,
            time_min=lambda wildcards, attempt: attempt * 20,
            tmpdir="tmp"
        log:
            "logs/mutect2/correct_fields/{sample}.log"
        params:
            rename_ad="'s/=AD;/=ADM;/g'",
            rename_ad_format="'s/:AD:/:ADM:/g'",
            fix_as_filterstatus="'s/ID=AS_FilterStatus,Number=A/ID=AS_FilterStatus,Number=1/g'"
        shell:
            "(gunzip -c {input} | "
            "sed {params.rename_ad} | "
            "sed {params.rename_ad_format} | "
            "sed {params.fix_as_filterstatus}) "
            "> {output} 2> {log}"

    ###############################
    ### Variant calling Mutect2 ###
    ###############################

    gatk_mutect2_somatic_config = {
        "genome": config["ref"]["fasta"],
        "known": config["ref"]["af_only"],
        "bed": config["ref"]["capture_kit_bed"],
        "dbsnp": config["ref"]["dbsnp"],
        "sample_list": design["Sample_id"].to_list()
    }


    module gatk_mutect2_somatic_meta:
        snakefile: "../../meta/bio/mutect2_somatic/test/Snakefile"
        config: gatk_mutect2_somatic_config


    use rule * from gatk_mutect2_somatic_meta




    ################################
    ### Variant Calling Varscan2 ###
    ################################

    varscan2_somatic_config = {
        "genome": config["ref"]["fasta"],
        "bed": config["ref"]["capture_kit_bed"]
    }

    module varscan2_somatic_meta:
        snakefile: "../../meta/bio/varscan2_somatic/test/Snakefile"
        config: varscan2_somatic_config

    use rule * from varscan2_somatic_meta


    ##############################
    ### GATK BAM RECALIBRATION ###
    ##############################

    gatk_bqsr_config = {
        "threads": config["threads"],
        "genome": config["ref"]["fasta"],
        "dbsnp": config["ref"]["dbsnp"]
    }

    module gatk_bqsr_meta:
        snakefile: "../../meta/bio/gatk_bqsr/test/Snakefile"
        config: gatk_bqsr_config


    use rule gatk_apply_baserecalibrator from gatk_bqsr_meta with:
        input:
            bam="picard/markduplicates/{sample}_{status}.bam",
            bam_index=get_bai("picard/markduplicates/{sample}_{status}.bam"),
            ref=config['ref']['fasta'],
            ref_idx=get_fai(config['ref']['fasta']),
            ref_dict=get_dict(config['ref']['fasta']),
            recal_table="gatk/recal_data_table/{sample}_{status}.grp"
        output:
            bam="gatk/recal_bam/{sample}_{status}.bam"
        message:
            "Applying BQSR on {wildcards.status} {wildcards.sample} with GATK"
        log:
            "logs/gatk/applybqsr/{sample}.{status}.log"


    use rule gatk_compute_baserecalibration_table from gatk_bqsr_meta with:
        input:
            bam="picard/markduplicates/{sample}_{status}.bam",
            bam_index=get_bai("picard/markduplicates/{sample}_{status}.bam"),
            ref=config['ref']['fasta'],
            ref_idx=get_fai(config['ref']['fasta']),
            ref_dict=get_dict(config['ref']['fasta']),
            known=config['ref']['dbsnp'],
            known_idx=get_tbi(config['ref']['dbsnp'])
        output:
            recal_table=temp("gatk/recal_data_table/{sample}_{status}.grp")
        message:
            "Compute BQSR table from {wildcards.status} {wildcards.sample} "
            "with GATK"
        log:
            "logs/gatk3/compute_bqsr/{sample}.{status}.log"


    #####################
    ### Deduplicating ###
    #####################

    rule picard_markduplicates:
        input:
            bam="samtools/sort/{sample}_{status}.bam",
            bai=get_bai("samtools/sort/{sample}_{status}.bam")
        output:
            bam=temp("picard/markduplicates/{sample}_{status}.bam"),
            metrics=temp(
                "picard/metrics/{sample}_{status}.picard.metrics.txt"
            )
        message:
            "Removing duplicates on {wildcards.sample} ({wildcards.status})"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 10240,
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir="tmp"
        log:
            "logs/picard/markduplicates/{sample}_{status}.log"
        params:
            "--ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES true"
        wrapper:
            "bio/picard/markduplicates"


    ###################
    ### BWA MAPPING ###
    ###################

    module bwa_fixmate_meta:
        snakefile: "../../meta/bio/bwa_fixmate/test/Snakefile"
        config: {"threads": config["threads"], "genome": config["ref"]["fasta"]}


    use rule samtools_index from bwa_fixmate_meta with:
        input:
            "samtools/sort/{sample}_{status}.bam"
        output:
            temp("samtools/sort/{sample}_{status}.bam.bai")
        message:
            "Indexing mapped reads of {wildcards.status} {wildcards.sample}"
        log:
            "logs/samtools/sort/{sample}.{status}.log"


    use rule samtools_sort_coordinate from bwa_fixmate_meta with:
        input:
            "samtools/fixmate/{sample}_{status}.bam"
        output:
            temp("samtools/sort/{sample}_{status}.bam")
        message:
            "Sorting {wildcards.status} {wildcards.sample} reads by query "
            "name for fixing mates"
        log:
            "logs/samtools/query_sort_{sample}.{status}.log"


    use rule samtools_fixmate from bwa_fixmate_meta with:
        input:
            "bwa_mem2/mem/{sample}_{status}.bam"
        output:
            temp("samtools/fixmate/{sample}_{status}.bam")
        message:
            "Fixing mate annotation on {wildcards.status} "
            "{wildcards.sample} with Samtools"
        log:
            "logs/samtools/fixmate/{sample}.{status}.log"


    use rule bwa_mem from bwa_fixmate_meta with:
        input:
            reads=expand(
                "fastp/trimmed/pe/{sample}_{status}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
            index=multiext(
                "bwa_mem2/index/genome", ".0123", ".amb", ".ann", ".pac"
            )
        output:
            temp("bwa_mem2/mem/{sample}_{status}.bam")
        message:
            "Mapping {wildcards.status} {wildcards.sample} with BWA"
        log:
            "logs/bwa_mem2/mem/{sample}.{status}.log"


    use rule bwa_index from bwa_fixmate_meta with:
        input:
            config["ref"]["fasta"]


    ############################
    ### FASTP FASTQ CLEANING ###
    ############################

    rule fastp_clean:
        input:
            sample=expand(
                "reads/{status}/{sample}.{stream}.fq.gz",
                stream=["1", "2"],
                allow_missing=True
            ),
        output:
            trimmed=expand(
                "fastp/trimmed/pe/{sample}_{status}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True
            ),
            html="fastp/html/pe/{sample}_{status}.fastp.html",
            json=temp("fastp/json/pe/{sample}_{status}.fastp.json")
        message: "Cleaning {wildcards.status} {wildcards.sample} with Fastp"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
            time_min=lambda wildcard, attempt: attempt * 45,
            tmpdir="tmp"
        params:
            adapters=config.get("fastp_adapters", None),
            extra=config.get("fastp_extra", "")
        log:
            "logs/fastp/{sample}.{status}.log"
        wrapper:
            "bio/fastp"


    #################################################
    ### Gather files from iRODS or mounting point ###
    #################################################

    rule bigr_copy:
        output:
            "reads/{status}/{sample}.{stream}.fq.gz"
        message:
            "Gathering {wildcards.status} {wildcards.sample} fastq files "
            "({wildcards.stream})"
        threads: 1
        resources:
            mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
            time_min=lambda wildcard, attempt: attempt * 45,
            tmpdir="tmp"
        params:
            input=lambda w, output: fastq_links[w.status][output[0]]
        log:
            "logs/bigr_copy/{status}/{sample}.{stream}.log"
        wrapper:
            "bio/BiGR/copy"


    ###########################
    ### Datasets indexation ###
    ###########################

    index_datasets_config = {
        "genome": config["ref"]["fasta"]
    }

    module index_datasets:
        snakefile: "../../meta/bio/index_datasets/test/Snakefile"
        config: index_datasets_config

    use rule * from index_datasets

    use rule samtools_bam_index from index_datasets with:
        input:
            "{command}/{subsommand}/{sample}_{status}.bam"
        output:
            get_bai("{command}/{subsommand}/{sample}_{status}.bam")
        message:
            "Indexing bam from {wildcards.command} ({wildcards.subcommand}) "
            "with samtools for {wildcards.sample} ({wildcards.status})"
        log:
            "logs/samtools/sort/{command}_{subsommand}/{sample}_{status}.log"




Authors
-------


* Thibault Dayris

* M boyba Diop

* Marc Deloger
