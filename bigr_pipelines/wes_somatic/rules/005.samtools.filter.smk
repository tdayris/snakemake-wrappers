rule samtools_filter_bed:
    input:
        "sambamba/sort/{sample}_{status}.bam",
        fasta=config["reference"]["fasta"],
        fasta_idx=get_fai(config["reference"]["fasta"]),
        fasta_dict=get_dict(config["reference"]["fasta"]),
        bed=config["reference"]["capture_kit_bed"],
    output:
        temp("samtools/filter/{sample}_{status}.bam"),
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    retries: 1
    params:
        extra=config["samtools"].get("view_filter", "-h -q 5"),
    log:
        "logs/samtools/filter/{sample}_{status}.log",
    wrapper:
        "bio/samtools/view"
