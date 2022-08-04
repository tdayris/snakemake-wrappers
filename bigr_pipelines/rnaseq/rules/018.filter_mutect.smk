rule muterc2_filter:
    input:
        vcf="mutect2/call/{sample}.vcf.gz",
        ref=config["reference"]["genome"],
        fasta_index=config["reference"]["genome_index"],
        fasta_dict=config["reference"]["genome_dict"],
        contamination="summary/{sample}_calculate_contamination.table",
        bam="gatk/splitncigarreads/{sample}.bam",
        bam_index="gatk/splitncigarreads/{sample}.bam.bai",
        f1r2="gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz",
    output:
        vcf=temp("mutect2/filter/{sample}.vcf.gz"),
    threads: 1
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["gatk"].get("filtermutectcalls", ""),
    log:
        "logs/mutect2/filter/{sample}.log",
    wrapper:
        "bio/gatk/filtermutectcalls"
