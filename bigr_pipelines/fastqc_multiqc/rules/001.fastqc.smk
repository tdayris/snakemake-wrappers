rule fastqc:
    input:
        "reads/{sample}.{stream}.fq.gz"
    output:
        html=temp("fastqc/{sample}.{stream}.html"),
        zip=temp("fastqc/{sample}_{stream}_fastqc.zip")
    message:
        "Assessing quality of {wildcards.sample}, ({wildcards.stream})"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024, 4096),
        time_min=lambda wildcard, attempt: attempt * 50,
        tmpdir="tmp"
    params:
        ""
    log:
        "logs/fastqc/{sample}.{stream}.log"
    wrapper:
        "bio/fastqc"