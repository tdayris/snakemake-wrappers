"""
This rule copies fastq files if needed, or symlinks
it also concatenates fastq files separated by commas in design
it also renames fastq files based on the sample name in design

No download is available now.
"""
rule bigr_copy:
    output:
        temp("data_input/{sample}.{stream}.fq.gz"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
    retries: 1
    params:
        input=lambda wildcards, output: fastq_links[output[0]],
    log:
        "logs/bigr_copy/{sample}.{stream}.log",
    wrapper:
        "bio/BiGR/copy"


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
        html=temp("fastp/html/{sample}.fastp.html"),
        json=temp("fastp/json/{sample}.fastp.json"),
    threads: min(config.get("max_threads", 5), 5)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        adapters=config["params"].get("fastp_adapters", None),
        extra=config["params"].get("fastp_extra", ""),
    log:
        "logs/fastp/{sample}.log",
    wrapper:
        "bio/fastp"
