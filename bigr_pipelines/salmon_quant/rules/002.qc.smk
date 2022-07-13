##############################################
### Clean and check quality of Fastq files ###
##############################################


rule fastp_clean:
    input:
        sample=expand(
            "reads/{sample}.{stream}.fq.gz", stream=["1", "2"], allow_missing=True
        ),
    output:
        trimmed=temp(
            expand(
                "fastp/trimmed/pe/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True,
            )
        ),
        html="fastp/html/pe/{sample}.fastp.html",
        json="fastp/json/pe/{sample}.fastp.json",
    message:
        "Cleaning {wildcards.sample} with Fastp"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
        time_min=lambda wildcard, attempt: attempt * 45,
        tmpdir="tmp",
    params:
        adapters=config["params"].get("fastp_adapters", None),
        extra=config["params"].get("fastp_extra", ""),
    log:
        "logs/fastp/{sample}.log",
    wrapper:
        "bio/fastp"


##########################################
### Check and assess origin of samples ###
##########################################


rule fastq_screen:
    input:
        "reads/{sample}.{stream}.fq.gz",
    output:
        txt=temp("fastq_screen/{sample}.{stream}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.{stream}.fastq_screen.png"),
    message:
        "Assessing quality of {wildcards.sample}, stream {wildcards.stream}"
    threads: config.get("threads", 20)
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024 * 8, 8192),
        time_min=lambda wildcard, attempt: attempt * 75,
        tmpdir="tmp",
    params:
        fastq_screen_config=config["fastq_screen"],
        subset=100000,
        aligner="bowtie2",
    log:
        "logs/fastq_screen/{sample}.{stream}.log",
    wrapper:
        "bio/fastq_screen"


#####################################################
### Gather and compiles multiple quality controls ###
#####################################################


rule multiqc_salmon_quant:
    input:
        salmon=expand(
            "salmon/pseudo_mapping/{sample}/quant.sf", sample=design["Sample_id"]
        ),
        html=expand("fastp/html/pe/{sample}.fastp.html", sample=design["Sample_id"]),
        json=expand("fastp/json/pe/{sample}.fastp.json", sample=design["Sample_id"]),
        fastq_screen=expand(
            "fastq_screen/{sample}.{stream}.fastq_screen.{ext}",
            sample=design["Sample_id"],
            stream=["1", "2"],
            ext=["txt", "png"],
        ),
    output:
        report(
            "multiqc/MultiQC.html",
            caption="../common/reports/multiqc.rst",
            category="Quality Controls",
        ),
    message:
        "Aggregating quality reports from Fastp and Salmon"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt * 1536, 10240),
        time_min=lambda wildcards, attempt: attempt * 35,
        tmpdir="tmp",
    log:
        "logs/multiqc.log",
    wrapper:
        "bio/multiqc"
