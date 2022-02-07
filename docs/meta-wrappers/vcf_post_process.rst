.. _`vcf_post_process`:

VCF_POST_PROCESS
================

Format gzipped VCF to TSV with renamed columns and human-friendly content when possible.


Example
-------

This meta-wrapper can be used by integrating the following into your workflow:

.. code-block:: python

    import sys
    from pathlib import Path

    worflow_source_dir = Path(next(iter(workflow.get_sources()))).absolute().parent
    common = str(worflow_source_dir / "../../../../bigr_pipelines/common/python")
    sys.path.append(common)

    from file_manager import *

    default_config_vcf_post_process = {
        "ncbi_build": "GRCh38",
        "center": "GustaveRoussy",
        "annotation_tag": "ANN=",
        "sample_list": list(),
        "genome": "/path/to/ref.fasta",
        "known": "/path/to/dbsnp",
        "gatk_filters": {},
        "chr": list(range(1, 23)) + ["X", "Y"]
    }

    try:
        if config == dict():
            config = default_config_vcf_post_process
    except NameError:
        config = default_config_vcf_post_process

    """
    Compress and index final VCF files
    """
    rule gath_final_vcf:
        input:
            expand(
                "maf/occurence_annotated/{sample}.vcf.gz",
                sample=config["sample_list"]
            ),
            expand(
                "maf/occurence_annotated/{sample}.vcf.gz.tbi",
                sample=config["sample_list"]
            )
        output:
            "final.vcf.list"
        message:
            "Aquiring list of final VCF files"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 128,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/somatic/list.log"
        shell:
            "echo {input} | sed 's/\s\+/\\n/g' > {output} 2> {log}"


    """
    Big concatenation for maftools
    """
    rule concat_all_mafs:
        input:
            expand(
                "maf/maftools/{sample}.maf",
                sample=config["sample_list"]
            )
        output:
            "maf/complete.maf"
        message:
            "Merging separates maf files in a single cohort maf"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir="tmp"
        log:
            "logs/maftools/concat.log"
        params:
            sed = "'s/FORMAT_{}_tumor_AF/Mutect2_Allele_Frequency/g'".format(
                config["sample_list"]
            )
        shell:
            "head -n 1 {input[0]} | sed {params.sed} > {output} 2> {log} && "
            "for VCF in {input}; do sed '1d' ${{VCF}}; done "
            ">> {output} 2>> {log}"


    """
    Assure MAFtools compatibility and human readability
    """
    rule rename_snpsift_maf_cols:
        input:
            tsv="maf/extracted/{sample}.tsv"
        output:
            tsv="maf/maftools/{sample}.maf"
        message:
            "Renaming columns to fit MAFtools requirement in {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024 * 5,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        group:
            "vcf_to_maf"
        log:
            "logs/maftools/rename/{sample}.log"
        params:
            add_cols=True,
            ncbi_build=config.get("NBCI_build", "GRCh38"),
            center=config.get("center", "GustaveRoussy"),
            Tumor_Sample_Barcode=lambda wildcards: f"{wildcards.sample}_tumor",
            Matched_Norm_Sample_Barcode=lambda wildcards: f"{wildcards.sample}_normal"
        wrapper:
            "bio/BiGR/rename_snpsift_maf_cols"




    """
    Extracting all INFO/FORMAT data, the list is built from vcf header
    """
    rule extract_all_fields:
        input:
            call="maf/occurence_annotated/{sample}.vcf"
        output:
            tsv=temp("maf/extracted/{sample}.tsv")
        message:
            "Extracting variant annotations for {wildcards.sample}"
        threads: 2
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 10240,
            time_min=lambda wildcards, attempt: attempt * 25,
            tmpdir="tmp"
        group:
            "vcf_to_maf"
        log:
            "logs/snpsift/extract_all_fields/{sample}.log"
        params:
            annotation_tag=config.get("annotation_tag", "ANN="),
            ignore_format=True,
            extra=config.get("extract_all_fields_extra", "-e '.' -s ','")
        wrapper:
            "bio/snpsift/extractAllFields"


    """
    Count variant occurence
    """
    rule variant_occurence_annotate:
        input:
            calls = ["maf/canonical/{sample}.vcf"],
            occurence = "maf/occurences.txt"
        output:
            calls = ["maf/occurence_annotated/{sample}.vcf"]
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/variant_occurence/uncompress/{sample}.log"
        wrapper:
            "bio/variantoccurence/annotate"


    rule concatenate_per_chr_information:
        input:
            expand("maf/{chr}/occurence.txt", chr=config["chr"])
        output:
            "maf/occurences.txt"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        log:
            "logs/variant_occurence/all.log"
        shell:
            "for i in {input}; do sed '1d' ${{i}}; done > {output} 2> {log}"


    rule variant_occurence_per_chr:
        input:
            calls=expand(
                "maf/canonical/{sample}.vcf.gz",
                sample=config["sample_list"]
            )
        output:
            txt="maf/{chr}/occurence.txt"
        threads: 7
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir="tmp"
        log:
            "logs/variant_occurence/{chr}.log"
        wrapper:
            "bio/variantoccurence/chromosomes"


    """
    Remove non-canonical chromosomes, and empty info fields
    """
    rule fix_vcf:
        input:
            vcf="maf/splitted/{sample}.vcf"
        output:
            vcf=temp("maf/canonical/{sample}.vcf")
        message:
            "Cleaning {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        group:
            "GLeaves"
        log:
            "logs/fix_vcf/{sample}.log"
        params:
            default_chr=config["chr"],
            remove_non_conventional_chromosomes=False
        wrapper:
            "bio/BiGR/fix_vcf"


    """
    Split annotation since it may lead to errors in MAFtools and/or in a result
    file more difficult to read by a human.
    """
    rule split_vcf_features:
        input:
            call="snpsift/format2info/{sample}.vcf"
        output:
            call=temp("maf/splitted/{sample}.vcf")
        message:
            "Splitting variant annotations for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        group:
            "GLeaves"
        log:
            "logs/split_vcf_features/{sample}.log"
        params:
            annotation_tag="ANN="
        wrapper:
            "bio/BiGR/split_vcf_features"


    """
    Copy format information, this is for end-users reading
    """
    rule format_to_info:
        input:
            call = "gatk/variant_filtration/{sample}.vcf"
        output:
            call = temp("snpsift/format2info/{sample}.vcf")
        message:
            "Moving format fields to info for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2048,
            time_min=lambda wildcards, attempt: attempt * 45,
            tmpdir="tmp"
        group:
            "GLeaves"
        log:
            "logs/vcf_format_to_info/{sample}.log"
        wrapper:
            "bio/BiGR/vcf_format_to_info"


    rule unzip_variant_filtration:
        input:
            "gatk/variant_filtration/{sample}.vcf.gz"
        output:
            temp("gatk/variant_filtration/{sample}.vcf")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        params:
            "-c"
        group:
            "GLeaves"
        log:
            "logs/gunzip/variant_filtration/{sample}"
        shell:
            "gunzip {params} {input} > {output} 2> {log}"


    """
    Variant calling quality control
    """
    rule gatk_variant_evaluation:
        input:
            vcf="gatk/variant_filtration/{sample}.vcf.gz",
            vcf_tbi=get_tbi("gatk/variant_filtration/{sample}.vcf.gz"),
            bam="picard/markduplicates/{sample}_tumor.bam",
            bai=get_bai("picard/markduplicates/{sample}_tumor.bam"),
            ref=config["genome"],
            fai=get_fai(config["genome"]),
            dict=get_dict(config["genome"]),
            known=config["known"],
            known_tbi=get_tbi(config["known"])
        output:
            directory("gatk_variant_evaluation/{sample}")
        message:
            "Evaluating variant calling of {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1024,
            time_min=lambda wildcards, attempt: attempt * 15,
            tmpdir="tmp"
        group:
            "GATK_Stats"
        log:
            "logs/gatk/varianteval/{sample}.log"
        params:
            extra=""
        wrapper:
            "bio/gatk/varianteval"


    """
    Add filter tags, these filters do not remove variants, only annotates
    """
    rule gatk_variant_filtration:
        input:
            vcf="snpsift/annotate_corrected/{sample}.vcf.gz",
            vcf_tbi=get_tbi("snpsift/annotate_corrected/{sample}.vcf.gz"),
            ref=config["genome"]
        output:
            vcf=temp("gatk/variant_filtration/{sample}.vcf.gz")
        message:
            "Filtering VCF for {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 10240,
            time_min=lambda wildcards, attempt: attempt * 25,
            tmpdir="tmp"
        group:
            "GATK_Stats"
        log:
            "logs/gatk/variant_filtration/{sample}.log"
        params:
            filters=config.get("gatk_filters", {
                "DepthBelow60X": "DP < 59",
                "BelowQualByDepth": "QD <= 2.0",
                "BelowBaseQuality": "QUAL < 30.0",
                "AboveFisherStrandBias": "FS > 60.0",
                "AboveStrandOddsRatio": "SOR > 3.0",
                "BelowMappingQuality": "MQ < 35.0",
                "BelowMQRankSum": "MQRankSum < -12.5",
                "BelowReadPosRankSum": "ReadPosRankSum < -8.0"
            })
        wrapper:
            "bio/gatk/variantfiltration"


    rule fix_annotation_for_gatk:
        input:
            call="snpsift/clinvar/{sample}.vcf"
        output:
            call=temp("snpsift/annotate_corrected/{sample}.vcf")
        message:
            "Correcting annotation type error in {wildcards.sample}"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 512,
            time_min=lambda wildcards, attempt: attempt * 25,
            tmpdir="tmp"
        params:
            remove_list=[
                "END=\([0-9]\+,\?\)\+"
            ],
            replace_dict=lambda wildcards: {
                ";;": ";",
                ";\t": "\t",
                ":ADM:": ":AD:",
                "ID=ADM,": "ID=AD,",
                f"FORMAT_{wildcards.sample}_tumor_AF": "Allele_Frequency"
            }
        group:
            "GATK_Stats"
        log:
            "logs/snpsift/annotate_corrected/{sample}.log"
        wrapper:
            "bio/sed"

Note that input, output and log file paths can be chosen freely, as long as the dependencies between the rules remain as listed here.
For additional parameters in each individual wrapper, please refer to their corresponding documentation (see links below).

When running with

.. code-block:: bash

    snakemake --use-conda

the software dependencies will be automatically deployed into an isolated environment before execution.



Used wrappers
---------------------

The following individual wrappers are used in this meta-wrapper:


* :ref:`bio/BiGR/split_vcf_features`

* :ref:`bio/snpsift/extractAllFields`

* :ref:`bio/BiGR/rename_snpsift_maf_cols`


Please refer to each wrapper in above list for additional configuration parameters and information about the executed code.






Notes
-----

VCF header will be used to extract information, no need to list them all.




Authors
-------


* Thibault Dayris

