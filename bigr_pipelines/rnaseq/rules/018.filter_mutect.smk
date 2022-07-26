rule muterc2_filter:
    input:
        vcf="mutect2/call/{sample}.vcf.gz",
        ref=config["genome"],
        fasta_index=get_fai(config["genome"]),
        fasta_dict=get_dict(config["genome"]),
        contamination="summary/{sample}_calculate_contamination.table",
        bam="sambamba/sort/{sample}.bam",
        bam_index=get_bai("sambamba/sort/{sample}.bam"),
        f1r2="gatk/artifacts_prior/{sample}.artifacts_prior.tar.gz",
    output:
        vcf=temp("mutect2/filter/{sample}.vcf.gz"),
    message:
        "Filtering Mutect2 calls for {wildcards.sample}"
    threads: 1
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["gatk"].get("filtermutectcalls", ""),
    log:
        "logs/mutect2/filter/{sample}.log",
    wrapper:
        "bio/gatk/filtermutectcalls"
