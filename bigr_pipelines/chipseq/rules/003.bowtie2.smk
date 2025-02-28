rule bowtie2_map:
    input:
        sample=expand(
            "fastp/trimmed/{sample}.{stream}.fastq", stream=streams, allow_missing=True
        ),
        idx=config["bowtie2"]["index"],
    output:
        temp("bowtie2/raw/{sample}.bam"),
    threads: 20
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_4h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bowtie2/{sample}.log",
    params:
        extra=(
            "--end-to-end "
            "--no-unal "
            "--no-mixed "
            "--no-discordant "
            "--phred33 "
            "-I 10 -X 700"
        ),
        index=f'{config["bowtie2"]["index_prefix"]}',
    wrapper:
        "bio/bowtie2/align"


rule sambamba_sort:
    input:
        mapping="bowtie2/raw/{sample}.bam",
    output:
        mapping=temp("bowtie2/sorted/{sample}.bam"),
    threads: 5
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    shadow:
        "shallow"
    log:
        "logs/sambamba/sort/{sample}.bwa.log",
    params:
        "",
    wrapper:
        "bio/sambamba/sort"


rule sambamba_index_bwa:
    input:
        "bowtie2/sorted/{sample}.bam",
    output:
        temp("bowtie2/sorted/{sample}.bam.bai"),
    threads: 8
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/sambamba/{sample}.index.bwa.log",
    params:
        extra="",
    wrapper:
        "bio/sambamba/index"
