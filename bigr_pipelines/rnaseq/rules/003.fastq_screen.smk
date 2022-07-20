"""
This snakefile calls fastqc on raw fastq files
"""

# Assess sample origin based on a wide range of potiential target genomes
rule fastq_screen:
    input:
        "data_input/{sample}.{stream}.fq.gz",
    output:
        txt=temp("fastq_screen/{sample}.{stream}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.{stream}.fastq_screen.png"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_75min_per_attempt,
        tmpdir="tmp",
    retries: 2
    params:
        fastq_screen_config=config["fastq_screen"]["genomes"],
        subset=config["fastq_screen"].get("subset", 100000),
        aligner=config["fastq_screen"].get("aligner", "bowtie2"),
    log:
        "logs/fastq_screen/{sample}.{stream}.log",
    wrapper:
        "bio/fastq_screen"