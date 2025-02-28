# Filter calls from mutect with estimated artifacts
"""
018.mutect2_filter
from
-> 017.learn_read_orientation_model
-> 017.calculate_contamination
-> 016.mutect2_germline
-> 010.gatk_split_n_cigar_reads
by
-> End job
"""


rule mutect2_filter:
    input:
        vcf="016.mutect2/call/{sample}.vcf.gz",
        ref=config["reference"]["genome"],
        fasta_index=config["reference"]["genome_index"],
        fasta_dict=config["reference"]["genome_dict"],
        contamination="017.summary/{sample}_calculate_contamination.table",
        bam="010.gatk/splitncigarreads/{sample}.bam",
        bam_index="010.gatk/splitncigarreads/{sample}.bam.bai",
        f1r2="010.gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz",
    output:
        vcf=temp("016.mutect2/filter/{sample}.vcf.gz"),
    threads: 1
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["gatk"].get("filtermutectcalls", ""),
    log:
        "logs/018.mutect2/filter/{sample}.log",
    wrapper:
        "bio/gatk/filtermutectcalls"
