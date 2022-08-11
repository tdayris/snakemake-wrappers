rule samtools_stats_raw:
    input:
        aln="bwa_mem2/sorted/{sample}_{status}.bam",
        bai="bwa_mem2/sorted/{sample}_{status}.bam.bai",
        ref=config["reference"]["fasta"],
        ref_idx=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        temp("samtools/stats/{sample}_{status}.raw.stats"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/samtools/stats/{sample}.{status}.raw.log",
    params:
        extra=config["samtools"].get("stats", ""),
    wrapper:
        "bio/samtools/stats"


rule samtools_stats_cleaned:
    input:
        aln="sambamba/markdup/{sample}_{status}.bam",
        bai="sambamba/markdup/{sample}_{status}.bam.bai",
        ref=config["reference"]["fasta"],
        ref_idx=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        temp("samtools/stats/{sample}_{status}.cleaned.stats"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/samtools/stats/{sample}.{status}.cleaned.log",
    params:
        extra=config["samtools"].get("stats", ""),
    wrapper:
        "bio/samtools/stats"


rule samtools_stats_fusions:
    input:
        aln="star/chimera/{sample}_{status}.bam",
        bai="star/chimera/{sample}_{status}.bam.bai",
        ref=config["reference"]["fasta"],
        ref_idx=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        temp("samtools/stats/{sample}_{status}.chimera.stats"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/samtools/stats/{sample}.{status}.chimera.log",
    params:
        extra=config["samtools"].get("stats", ""),
    wrapper:
        "bio/samtools/stats"


rule fastq_screen:
    input:
        "fastp/trimmed/{sample}_{status}.{stream}.fastq",
    output:
        txt=temp("fastq_screen/{sample}.{stream}.{status}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.{stream}.{status}.fastq_screen.png"),
    threads: config.get("threads", 20)
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_75min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        fastq_screen_config=config["fastq_screen"],
        subset=100000,
        aligner="bowtie2",
    log:
        "logs/fastqc/{sample}.{stream}.{status}.log",
    wrapper:
        "bio/fastq_screen"
