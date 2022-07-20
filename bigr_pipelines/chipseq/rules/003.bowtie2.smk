rule bowtie2_map:
    input:
        sample=expand(
            "fastp/trimmed/{sample}.{stream}.fastq", sample=sample_list, stream=streams
        ),
        idx=config["bowtie2"]["index"],
    output:
        temp("bowtie2/raw/{sample}.bam"),
    threads: 20
    resources:
        mem_mb=get_75gb_and_5gb_per_attempt,
        time_min=get_4h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bowtie2/{sample}.log",
    params:
        extra="--end-to-end --very-sensitive --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700",
    wrapper:
        "bio/bowtie2/align"


rule samtools_index_bwa:
    input:
        "bowtie2/raw/{sample}.bam",
    output:
        temp("bowtie2/raw/{sample}.bam.bai"),
    threads: 4
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "log/samtools/{sample}.index.bwa.log",
    params:
        extra="",
    wrapper:
        "bio/samtools/index"
