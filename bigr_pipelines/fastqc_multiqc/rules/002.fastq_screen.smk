rule fastq_screen:
    input:
        "reads/{sample}.{stream}.fq.gz"
    output:
        txt=temp("fastq_screen/{sample}.{stream}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.{stream}.fastq_screen.png")
    message:
        "Assessing quality of {wildcards.sample}, stream {wildcards.stream}"
    threads: config.get("threads", 20)
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 4096, 8192),
        time_min=lambda wildcard, attempt: attempt * 50,
        tmpdir="tmp"
    params:
        fastq_screen_config=config["fastq_screen"],
        subset=100000,
        aligner='bowtie2'
    log:
        "logs/fastq_screen/{sample}.{stream}.log"
    wrapper:
        str(worflow_source_dir /  "bio" / "fastq_screen")