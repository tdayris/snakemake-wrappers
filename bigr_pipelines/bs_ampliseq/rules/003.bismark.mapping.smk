rule bismark_mapping:
    input:
        genome=config["ref"]["fasta"],
        fq_1="fastp/trimmed/{sample}.1.fastq",
        fq_2="fastp/trimmed/{sample}.2.fastq",
        bismark_indexes_dir=config["ref"]["bismark_index"],
        genomic_freq=config["ref"]["bismark_frequencies"],
    output:
        bam=temp("bismark/align/{sample}.PE.bam"),
        report=temp("bismark/align/{sample}_PE_report.txt"),
        nucleotide_stats=temp("bismark/align/{sample}.nucleotide_stats.txt"),
        bam_unmapped_1=temp("bismark/align/{sample}_unmapped_reads_1.fq.gz"),
        bam_unmapped_2=temp("bismark/align/{sample}_unmapped_reads_2.fq.gz"),
        ambiguous_1=temp("bismark/align/{sample}_ambiguous_reads_1.fq.gz"),
        ambiguous_2=temp("bismark/align/{sample}_ambiguous_reads_2.fq.gz"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_75gb_and_2gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bismark/mapping/{sample}.log",
    params:
        extra=lambda wildcards: f"{config['bismark']['mapping_extra']} --rg_tag '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPU:{wildcards.sample}\\tPL:ILLUMINA\\tCN:IGR\\tDS:BSAmpliSeq\\tPG:BOWTIE2'",
    wrapper:
        "bio/bismark/bismark"
