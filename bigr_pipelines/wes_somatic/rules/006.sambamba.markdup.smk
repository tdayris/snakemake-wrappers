rule sambamba_markduplicates:
    input:
        bam="samtools/filter/{sample}_{status}.bam",
        bai="samtools/filter/{sample}_{status}.bam.bai",
    output:
        bam=temp("sambamba/markdup/{sample}_{status}.bam"),
    threads: 10
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "logs/sambamba/markduplicates/{sample}_{status}.log",
    params:
        extra=config["sambamba"].get("markdup_extra", "--remove-duplicates"),
    wrapper:
        "bio/sambamba/markdup"
