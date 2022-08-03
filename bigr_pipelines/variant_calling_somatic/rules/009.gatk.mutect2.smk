##################################
### Correct multiallelic sites ###
### and annotation issues      ###
##################################


rule correct_mutect2_vcf:
    input:
        "bcftools/mutect2/{sample}.vcf.gz",
    output:
        temp("mutect2/corrected/{sample}.vcf"),
    threads: 2
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_20min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/mutect2/correct_fields/{sample}.log",
    params:
        fix_as_filterstatus="'s/ID=AS_FilterStatus,Number=A/ID=AS_FilterStatus,Number=1/g'",
    shell:
        "(gunzip -c {input} | "
        "sed {params.fix_as_filterstatus}) "
        "> {output} 2> {log}"


rule split_multiallelic_mutect2:
    input:
        call="mutect2/filter/{sample}.vcf.gz",
        idx=get_tbi("mutect2/filter/{sample}.vcf.gz"),
        fasta=config["reference"]["fasta"],
        fasta_idx=config["reference"]["fasta_index"],
        fasta_dict=config["reference"]["fasta_dict"],
    output:
        "bcftools/mutect2/{sample}.vcf.gz",
    threads: 2
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["bcftools"].get("split_multiallelic", "-m -both --check-ref w"),
    log:
        "logs/bcftools/norm/mutect2/{sample}.log",
    wrapper:
        "bio/bcftools/norm"


###########################################
### Estimate cross-sample contamination ###
###########################################


rule gatk_filter_mutect_calls:
    input:
        vcf="mutect2/call/{sample}.vcf.gz",
        vcf_tbi=get_tbi("mutect2/call/{sample}.vcf.gz"),
        ref=config["reference"]["fasta"],
        ref_index=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
        bam="sambamba/markdup/{sample}_tumor.bam",
        bam_index=get_bai("sambamba/markdup/{sample}_tumor.bam"),
        f1r2="gatk/orientation_model/{sample}/{sample}.artifacts-prior.tar.gz",
        contamination="summary/{sample}_calculate_contamination.table",
    output:
        vcf=temp("mutect2/filter/{sample}.vcf.gz"),
        vcf_index=temp("mutect2/filter/{sample}.vcf.gz.tbi"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 2
    params:
        extra=config["gatk"].get(
            "filter_mutect_calls",
            "--create-output-variant-index --min-median-mapping-quality 35",
        ),
    log:
        "logs/mutect2/filter/{sample}.log",
    wrapper:
        "bio/gatk/filtermutectcalls"


"""
Estimate possible contaminations
"""


rule calculate_tumor_contamination:
    input:
        summary="gatk/getpileupsummaries/{sample}_tumor_getpileupsummaries.table",
        normal="gatk/getpileupsummaries/{sample}_normal_getpileupsummaries.table",
    output:
        table=temp("summary/{sample}_calculate_contamination.table"),
        segmentation=temp("summary/{sample}_segments.table"),
    group:
        "Contamination_Estimate"
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["gatk"].get("calculate_contamination", ""),
    log:
        "logs/gatk/CalculateContamination/{sample}.log",
    wrapper:
        "bio/gatk/calculatecontamination"


"""
Summarize the read support over known variants
"""


rule get_pileup_summaries:
    input:
        bam="sambamba/markdup/{sample}_{status}.bam",
        bam_index=get_bai("sambamba/markdup/{sample}_{status}.bam"),
        intervals=config["reference"]["capture_kit_bed"],
        variants=config["reference"]["af_only"],
        variants_index=config["reference"]["af_only_tbi"],
    output:
        table=temp("gatk/getpileupsummaries/{sample}_{status}_getpileupsummaries.table"),
    group:
        "Contamination_Estimate"
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["gatk"].get("get_pileup_summaries", ""),
    log:
        "logs/gatk/GetPileupSummaries/{sample}.{status}.log",
    wrapper:
        "bio/gatk/getpileupsummaries"


"""
Build orientation bias model to filter false positive calls
"""


rule gatk_learn_read_orientation_model:
    input:
        f1r2="mutect2/f1r2/{sample}.tar.gz",
    output:
        temp("gatk/orientation_model/{sample}/{sample}.artifacts-prior.tar.gz"),
    threads: 1
    resources:
        time_min=get_1h_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["gatk"].get("learn_read_orientation_model", ""),
    log:
        "gatk/orientation_model/{sample}.log",
    wrapper:
        "bio/gatk/learnreadorientationmodel"


######################
### Actual Calling ###
######################

"""
This rule calls somatic variants with GATK Mutect2
"""


rule mutect2_somatic:
    input:
        fasta=config["reference"]["fasta"],
        fasta_idx=config["reference"]["fasta_index"],
        fasta_dict=config["reference"]["fasta_dict"],
        map="sambamba/markdup/{sample}_normal.bam",
        map_idx="sambamba/markdup/{sample}_normal.bam.bai",
        tumor="sambamba/markdup/{sample}_tumor.bam",
        tumor_idx="sambamba/markdup/{sample}_tumor.bam.bai",
        germline=config["reference"]["af_only"],
        germline_tbi=config["reference"]["af_only_tbi"],
        intervals=config["reference"]["capture_kit_bed"],
    output:
        vcf=temp("mutect2/call/{sample}.vcf.gz"),
        vcf_index=temp("mutect2/call/{sample}.vcf.gz.tbi"),
        f1r2=temp("mutect2/f1r2/{sample}.tar.gz"),
    threads: config.get("max_threads", 20)
    resources:
        time_min=get_5h_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir="tmp",
    retries: 3
    params:
        extra=lambda wildcards, output: (
            f"{config['gatk'].get('mutect2', '')} "
            f"--tumor-sample {wildcards.sample}_tumor "
            f"--normal-sample {wildcards.sample}_normal "
        ),
    log:
        "logs/gatk/mutect2/call/{sample}.log",
    wrapper:
        "bio/gatk/mutect"
