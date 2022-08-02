################################
### Estimate sequencing bias ###
################################
"""
Build orientation model from f1r2 counts made in Mutect2
"""


rule learn_read_orientation_model:
    input:
        f1r2="mutect2/f1r2/{sample}.tar.gz",
    output:
        temp("gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["gatk"].get("learn_read_orientation_model", ""),
    log:
        "logs/gatk/learnreadorientationmodel/{sample}.log",
    wrapper:
        "bio/gatk/learnreadorientationmodel"


###########################################
### Estimate cross-sample contamination ###
###########################################


"""
Estimate possible contaminations
"""


rule calculate_contamination:
    input:
        summary="gatk/getpileupsummaries/{sample}_getpileupsummaries.table",
    output:
        table=temp("summary/{sample}_calculate_contamination.table"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
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
        bam="gatk/splitncigarreads/{sample}.bam",
        bam_index="gatk/splitncigarreads/{sample}.bam.bai",
        intervals=config["reference"]["capturekit_bed"],
        variants=config["reference"]["af_only"],
        variants_index=config["reference"]["af_only"],
    output:
        table=temp("gatk/getpileupsummaries/{sample}_getpileupsummaries.table"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["gatk"].get("pileup_summaries", ""),
    log:
        "logs/gatk/GetPileupSummaries/{sample}.log",
    wrapper:
        "bio/gatk/getpileupsummaries"
