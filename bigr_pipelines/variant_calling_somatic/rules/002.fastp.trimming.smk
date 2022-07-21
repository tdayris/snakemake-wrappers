rule fastp_clean:
    input:
        sample=expand(
            "reads/{status}/{sample}.{stream}.fq.gz",
            stream=["1", "2"],
            allow_missing=True,
        ),
    output:
        trimmed=temp(
            expand(
                "fastp/trimmed/pe/{sample}_{status}.{stream}.fastq",
                stream=["1", "2"],
                allow_missing=True,
            )
        ),
        html="fastp/html/pe/{sample}_{status}.fastp.html",
        json="fastp/json/pe/{sample}_{status}.fastp.json",
    message:
        "Cleaning {wildcards.status} {wildcards.sample} with Fastp"
    threads: 10
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 4096, 15360),
        time_min=lambda wildcard, attempt: attempt * 45,
        tmpdir="tmp",
    params:
        adapters=config.get("fastp_adapters", None),
        extra=config.get("fastp_extra", ""),
    log:
        "logs/fastp/{sample}.{status}.log",
    wrapper:
        "bio/fastp"
