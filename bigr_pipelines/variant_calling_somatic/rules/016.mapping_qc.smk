rule alignment_summary:
    input:
        bam="sambamba/sort/{sample}_{status}.bam",
        bam_index=get_bai("sambamba/sort/{sample}_{status}.bam"),
        ref=config['ref']['fasta'],
        ref_idx=get_fai(config['ref']['fasta']),
        ref_dict=get_dict(config['ref']['fasta']),
    output:
        temp("picard/alignment_summary/{sample}_{status}.summary.txt")
    message:
        "Collecting alignment metrics on GATK recalibrated {wildcards.sample}"
        " (considering {wildcards.status})"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020,
        time_min=lambda wildcards, attempt: attempt * 45,
        tmpdir="tmp"
    log:
        "logs/picard/alignment_summary/{sample}_{status}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    wrapper:
        "bio/picard/collectalignmentsummarymetrics"


rule fastq_screen:
    input:
        "reads/{status}/{sample}.{stream}.fq.gz"
    output:
        txt=temp("fastq_screen/{sample}.{stream}.{status}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.{stream}.{status}.fastq_screen.png")
    message:
        "Assessing quality of {wildcards.sample}, {wildcards.stream}"
        " (considering {wildcards.status})"
    threads: config.get("threads", 20)
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024 * 8, 20480),
        time_min=lambda wildcard, attempt: attempt * 75,
        tmpdir="tmp"
    params:
        fastq_screen_config=config["fastq_screen"],
        subset=100000,
        aligner='bowtie2'
    log:
        "logs/fastqc/{sample}.{stream}.{status}.log"
    wrapper:
        "bio/fastq_screen"