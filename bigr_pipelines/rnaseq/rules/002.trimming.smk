"""
This snakefile handles trimming and QC on raw fastq files
"""


# Clean and check quality of Fastq files
"""
002.fastp_clean
from
-> 001.bigr_copy
by
-> 003.fastq_screen
-> 004.salmon_quant
"""


rule fastp_clean:
    input:
        sample=expand(
            "data_input/{sample}.{stream}.fq.gz", stream=["1", "2"], allow_missing=True
        ),
    output:
        trimmed=temp(
            expand(
                "002.fastp/trimmed/{sample}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True,
            )
        ),
        html=temp("002.fastp/html/{sample}.fastp.html"),
        json=temp("002.fastp/json/{sample}.fastp.json"),
    threads: min(config.get("max_threads", 5), 5)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        adapters=config["fastp"].get("fastp_adapters", None),
        extra=config["fastp"].get("fastp_extra", ""),
    log:
        "logs/002.fastp/{sample}.log",
    wrapper:
        "bio/fastp"
