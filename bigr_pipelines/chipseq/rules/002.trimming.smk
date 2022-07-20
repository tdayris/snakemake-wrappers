"""
This snakefile handles trimming and QC on raw fastq files
"""


# Clean and check quality of Fastq files
rule fastp_clean:
    input:
        sample=expand(
            "data_input/{sample}.{stream}.fq.gz", stream=["1", "2"], allow_missing=True
        ),
    output:
        trimmed=temp(
            expand(
                "fastp/trimmed/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True,
            )
        ),
        html="fastp/html/{sample}.fastp.html",
        json="fastp/json/{sample}.fastp.json",
    threads: min(config.get("max_threads", 5), 5)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    params:
        adapters=config["fastp"].get("fastp_adapters", None),
        extra=config["fastp"].get("fastp_extra", ""),
    log:
        "logs/fastp/{sample}.log",
    wrapper:
        "bio/fastp"
