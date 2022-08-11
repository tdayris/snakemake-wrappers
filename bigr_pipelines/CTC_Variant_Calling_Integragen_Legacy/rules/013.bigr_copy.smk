"""
This Snakefile deals with IO from and to iRODS
"""


# Copy files on BiGR Flamingo
# iRODS paths are accepted
rule bigr_copy:
    output:
        temp("data_input/{sample}.{status}.bam"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
    retries: 1
    params:
        input=lambda wildcards, output: bam_links[output[0]],
    log:
        "logs/bigr_copy/{sample}.{status}.log",
    wrapper:
        "bio/BiGR/copy"
