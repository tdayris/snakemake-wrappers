# Trimm and perform MultiQC complient QC
rule fastp_clean:
    input:
        sample=expand(
            "reads/{sample}.{stream}.fq.gz",
            stream=["1", "2"],
            allow_missing=True
        ),
    output:
        trimmed=expand(
            "fastp/trimmed/pe/{sample}.{stream}.fastq",
            stream=["1", "2"],
            allow_missing=True
        ),
        html="fastp/html/pe/{sample}.fastp.html",
        json=temp("fastp/json/pe/{sample}.fastp.json")
    message: "Cleaning {wildcards.sample} with Fastp"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
        time_min=lambda wildcard, attempt: attempt * 45
    params:
        adapters=config.get("fastp_adapters", None),
        extra=config.get("fastp_extra", "")
    log:
        "logs/fastp/{sample}.log"
    wrapper:
        "bio/fastp"


# Gather all QC reports
rule multiqc_ampliseq:
    input:
        html=expand(
            "fastp/html/pe/{sample}.fastp.html",
            sample=design["Sample_id"]
        ),
        json=expand(
            "fastp/json/pe/{sample}.fastp.json",
            sample=design["Sample_id"]
        ),
        picard=expand(
            "picard/alignment_summary/{sample}.summary.txt",
            sample=design["Sample_id"]
        ),
        fastq_screen=expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
            sample=design["Sample_id"],
            stream=["1", "2"],
            ext=["txt", "png"]
        )
    output:
        report(
            "multiqc/variant_calling_ampliseq.html",
            caption="../common/reports/multiqc.rst",
            category="Quality Controls"
        )
    message:
        "Aggregating quality reports from SnpEff"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
        time_min=lambda wildcards, attempt: attempt * 35
    log:
        "logs/multiqc.log"
    wrapper:
        "bio/multiqc"


# Check bwa mapping
rule alignment_summary:
    input:
        bam="sambamba/sort/{sample}.bam",
        bam_index="sambamba/sort/{sample}.bam.bai",
        ref=config['ref']['fasta'],
        ref_idx=get_fai(config['ref']['fasta']),
        ref_dict=get_dict(config['ref']['fasta']),
    output:
        temp("picard/alignment_summary/{sample}.summary.txt")
    message:
        "Collecting alignment metrics on GATK recalibrated {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1020,
        time_min=lambda wildcards, attempt: attempt * 45
    log:
        "logs/picard/alignment_summary/{sample}.log"
    params:
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    wrapper:
        "bio/picard/collectalignmentsummarymetrics"


# Check organism origin and assess contaminations
rule fastq_screen:
    input:
        "reads/{sample}.{stream}.fq.gz"
    output:
        txt=temp("fastq_screen/{sample}.{stream}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.{stream}.fastq_screen.png")
    message:
        "Assessing quality of {wildcards.sample}, {wildcards.stream}"
    threads: config.get("threads", 20)
    resources:
        mem_mb=lambda wildcard, attempt: attempt * 1024 * 8,
        time_min=lambda wildcard, attempt: attempt * 50
    params:
        fastq_screen_config=config["fastq_screen"],
        subset=100000,
        aligner='bowtie2'
    log:
        "logs/fastqscreen/{sample}.{stream}.log"
    wrapper:
        "bio/fastq_screen"