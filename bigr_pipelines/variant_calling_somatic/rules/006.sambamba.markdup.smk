rule sambamba_markduplicates:
    input:
        bam="samtools/filter/{sample}_{status}.bam",
        bai=get_bai("samtools/filter/{sample}_{status}.bam")
    output:
        bam=temp("sambamba/markdup/{sample}_{status}.bam")
    message:
        "Removing duplicates on {wildcards.sample} ({wildcards.status})"
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10240,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/sambamba/markduplicates/{sample}_{status}.log"
    params:
        extra = config.get(
            "sambamba", {"markdup": "--remove-duplicates"}
        ).get("markdup", "--remove-duplicates")
    wrapper:
        "bio/sambamba/markdup"