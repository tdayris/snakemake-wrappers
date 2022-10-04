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
        "reads/{sample}.fq.gz"
    output:
        html=temp("fastqc/{sample}.html"),
        zip=temp("fastqc/{sample}_fastqc.zip")
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    params:
        ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "bio/fastqc"