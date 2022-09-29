# Contamination quality control
"""
001.fastq_screen
from
-> 000.bigr_copy
by
-> 003.multiqc
-> 003.irods_complient
"""
rule fastq_screen:
    input:
        "reads/{sample}.fq.gz"
    output:
        txt=temp("fastq_screen/{sample}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.fastq_screen.png")
    threads: config.get("threads", 20)
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024 * 10),
        time_min=lambda wildcard, attempt: attempt * 50,
        tmpdir="tmp"
    params:
        fastq_screen_config=config["fastq_screen"],
        subset=100000,
        aligner='bowtie2'
    log:
        "logs/fastq_screen/{sample}.log"
    wrapper:
        "bio/fastq_screen"