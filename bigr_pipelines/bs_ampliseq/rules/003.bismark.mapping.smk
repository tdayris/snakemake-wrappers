rule bismark_mapping:
    input:
        genome=config["ref"]["fasta"],
        fq_1="fastp/trimmed/{sample}.1.fastq",
        fq_2="fastp/trimmed/{sample}.2.fastq",
        bismark_indexes_dir=config["ref"]["bismark_index"],
        genomic_freq=config["ref"]["bismark_frequencies"],
    output:
        bam=temp("bismark/align/{sample}.bam"),
        report=temp("bismark/align/{sample}_PE_report.txt"),
        nucleotide_stats=temp("bismark/align/{sample}.nucleotide_stats.txt"),
    threads: 4
    resources:
        mem_mb=get_75gb_and_2gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bismark/mapping/{sample}.log",
    params:
        extra=lambda wildcards: f"{config['bismark']['mapping_extra']} --rg_id {wildcards.sample} --rg_sample {wildcards.sample}",
    wrapper:
        "bio/bismark/bismark"
