# Quality control
"""
001.fastqc
from
-> 000.bigr_copy
by
-> 003.multiqc
-> 003.irods_complient
"""
rule fastqc:
    input:
        "reads/{sample}.{stream}.fq.gz"
    output:
        html=temp("fastqc/{sample}.html"),
        zip=temp("fastqc/{sample}_fastqc.zip")
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