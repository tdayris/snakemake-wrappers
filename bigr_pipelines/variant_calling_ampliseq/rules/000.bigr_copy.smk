# Gather files from design paths
# iRODS paths in Flamingo are accepted
rule bigr_copy_fq:
    output:
        "reads/{sample}.{stream}.fq.gz"
    message:
        "Gathering {wildcards.sample} fastq file ({wildcards.stream})"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
        time_min=lambda wildcard, attempt: attempt * 45
    params:
        input=lambda wildcards, output: fastq_links[output[0]]
    log:
        "logs/bigr_copy/{sample}.{stream}.log"
    wrapper:
        "bio/BiGR/copy"