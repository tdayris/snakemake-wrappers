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
        "reads/{sample}.fq.gz",
    output:
        txt=temp("fastq_screen/{sample}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.fastq_screen.png"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_75min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        fastq_screen_config=config["fastq_screen"]["database"],
        subset=config["fastq_screen"].get("subset", 100000),
        aligner=config["fastq_screen"].get("aligner", "bowtie2"),
    log:
        "logs/003.fastq_screen/{sample}.log",
    wrapper:
        "bio/fastq_screen"