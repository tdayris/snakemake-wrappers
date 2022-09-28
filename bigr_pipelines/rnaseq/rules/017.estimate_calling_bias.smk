################################
### Estimate sequencing bias ###
################################
# Build orientation model from f1r2 counts made in Mutect2
"""
017.learn_read_orientation_model
from
-> 016.mutect2_germline
by
-> 018.mutect2_filter
"""


rule learn_read_orientation_model:
    input:
        f1r2="016.mutect2/f1r2/{sample}.tar.gz",
    output:
        temp("010.gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["gatk"].get("learn_read_orientation_model", ""),
    log:
        "logs/017.gatk/learnreadorientationmodel/{sample}.log",
    wrapper:
        "bio/gatk/learnreadorientationmodel"


###########################################
### Estimate cross-sample contamination ###
###########################################
# Estimate possible contaminations
"""
017.calculate_contamination
from
-> 017.calculate_contamination
by
-> 018.mutect2_filter
"""


rule calculate_contamination:
    input:
        summary="g010.atk/getpileupsummaries/{sample}_getpileupsummaries.table",
    output:
        table=temp("017.summary/{sample}_calculate_contamination.table"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["gatk"].get("calculate_contamination", ""),
    log:
        "logs/017.gatk/CalculateContamination/{sample}.log",
    wrapper:
        "bio/gatk/calculatecontamination"


# Summarize the read support over known variants
"""
017.get_pileup_summaries
from
-> 010.gatk_split_n_cigar_reads
by
-> 017.calculate_contamination
"""


rule get_pileup_summaries:
    input:
        bam="010.gatk/splitncigarreads/{sample}.bam",
        bam_index="010.gatk/splitncigarreads/{sample}.bam.bai",
        intervals=config["reference"]["capturekit_bed"],
        variants=config["reference"]["af_only"],
        variants_index=config["reference"]["af_only"],
    output:
        table=temp("010.gatk/getpileupsummaries/{sample}_getpileupsummaries.table"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["gatk"].get("pileup_summaries", ""),
    log:
        "logs/017.gatk/GetPileupSummaries/{sample}.log",
    wrapper:
        "bio/gatk/getpileupsummaries"
