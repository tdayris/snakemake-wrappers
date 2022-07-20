rule bedtools_genomecoveragebam:
    input:
        "samtools/view/{sample}.bam",
        "samtools/view/{sample}.bam.bai",
    output:
        temp("bedtools/genomecov/{sample}.bedgraph"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bedtools/genomecov/{sample}.log",
    params:
        "-bg",
    wrapper:
        "bio/bedtools/genomecov"
