#####################
### Apply filters ###
#####################

"""
005::mutect2_filter:
-> 005::learn_read_orientation_model
-> 005::calculate_contamination
-> 005::mutect2_germline
-> 004::samtools_view_filter
-> 004::sambamba_index
-> 001::samtools_index_genome
-> 001::picard_createsequencedictionary

Filter over estimated contaminations
"""


rule mutect2_filter:
    input:
        vcf="mutect2/call/{sample}.vcf.gz",
        ref=config[genome_id]["fasta"],
        fasta_index=fai_file,
        fasta_dict=dict_file,
        contamination="summary/{sample}_calculate_contamination.table",
        bam="samtools/filter/{sample}.bam",
        bam_idx="samtools/filter/{sample}.bam.bai",
        f1r2="gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz",
    output:
        vcf=temp("mutect2/filter/{sample}.vcf.gz"),
        vcf_tbi=temp("mutect2/filter/{sample}.vcf.gz.tbi"),
    threads: 1
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_10gb_per_attempt,
        tmpdir="tmp",
    params:
        extra=("--create-output-variant-index"),
    log:
        "logs/mutect2/filter/{sample}.log",
    wrapper:
        "bio/gatk/filtermutectcalls"


"""
005::learn_read_orientation_model:
-> 005::mutect2_germline

Build orientation model from f1r2 counts made in Mutect2
"""


rule learn_read_orientation_model:
    input:
        f1r2="mutect2/f1r2/{sample}.tar.gz",
    output:
        temp("gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    params:
        extra="",
    log:
        "logs/gatk/learnreadorientationmodel/{sample}.log",
    wrapper:
        "bio/gatk/learnreadorientationmodel"


"""
005::calculate_contamination:
-> 005::get_pileup_summaries

Estimate possible contaminations
"""


rule calculate_contamination:
    input:
        summary="gatk/getpileupsummaries/{sample}_getpileupsummaries.table",
    output:
        table=temp("summary/{sample}_calculate_contamination.table"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    params:
        extra="",
    log:
        "logs/gatk/CalculateContamination/{sample}.log",
    wrapper:
        "bio/gatk/calculatecontamination"


"""
005::get_pileup_summaries:
-> 004::samtools_view_filter
-> 004::sambamba_index
-> 001::tabix_index

Summarize the read support over known variants
"""


rule get_pileup_summaries:
    input:
        bam="samtools/filter/{sample}.bam",
        bam_idx="samtools/filter/{sample}.bam.bai",
        intervals=config[genome_id]["bed"],
        variants=config[genome_id]["vcf"],
        variants_index=tbi_file,
    output:
        table=temp("gatk/getpileupsummaries/{sample}_getpileupsummaries.table"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    params:
        extra="",
    log:
        "logs/gatk/GetPileupSummaries/{sample}.log",
    wrapper:
        "bio/gatk/getpileupsummaries"


"""
005::mutect2_germline:
-> 001::samtools_index_genome
-> 001::picard_createsequencedictionary
-> 004::samtools_view_filter
-> 004::sambamba_index
-> 001::tabix_index

This rule calls germline variants with GATK Mutect2

# https://gatk.broadinstitute.org/hc/en-us/articles/360042479112-CreateSomaticPanelOfNormals-BETA-
# Note that as of May, 2019 -max-mnp-distance 
# must be set to zero to avoid a bug in GenomicsDBImport.
"""


rule mutect2_germline:
    input:
        fasta=config[genome_id]["fasta"],
        fasta_index=fai_file,
        fasta_dict=dict_file,
        map="samtools/filter/{sample}.bam",
        bam_idx="samtools/filter/{sample}.bam.bai",
        germline=config[genome_id]["vcf"],
        germline_tbi=tbi_file,
        intervals=config[genome_id]["bed"],
    output:
        vcf=temp("mutect2/call/{sample}.vcf.gz"),
        f1r2=temp("mutect2/f1r2/{sample}.tar.gz"),
    threads: 10
    resources:
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
        mem_mb=get_15gb_per_attempt,
    params:
        extra=(
            "--max-mnp-distance 0 "
            "--max-reads-per-alignment-start 0 "
            "--tumor-sample Mutect2_{sample} "
            "--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter "
        ),
    log:
        "logs/gatk/mutect2/call/{sample}.log",
    wrapper:
        "bio/gatk/mutect"
