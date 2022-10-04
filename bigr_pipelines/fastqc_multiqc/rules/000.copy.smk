# Copies / Links files from irods or fs path
"""
000.bigr_copy
from
-> Entry job
by
-> 001.fastqc
-> 002.fastq_screen
"""
rule bigr_copy:
    output:
        temp("reads/{sample}.fq.gz")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    params:
        input=lambda wildcards, output: fastq_links[output[0]]
    log:
        "logs/bigr_copy/{sample}.log"
    wrapper:
        "bio/BiGR/copy"