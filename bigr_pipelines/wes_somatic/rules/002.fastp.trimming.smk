rule fastp_clean:
    input:
        sample=expand(
            "data_input/{status}/{sample}.{stream}.fq.gz",
            stream=["1", "2"],
            allow_missing=True,
        ),
    output:
        trimmed=temp(
            expand(
                "fastp/trimmed/{sample}_{status}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True,
            )
        ),
        html="fastp/html/pe/{sample}_{status}.fastp.html",
        json="fastp/json/pe/{sample}_{status}.fastp.json",
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        adapters=config["fastp"].get("fastp_adapters", None),
        extra=config["fastp"].get("fastp_extra", ""),
    log:
        "logs/fastp/{sample}.{status}.log",
    wrapper:
        "bio/fastp"
