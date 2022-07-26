rule arriba:
    input:
        bam="star/{sample}/{sample}.bam",
        genome=config["reference"]["genome"],
        annotation=config["reference"]["gtf"],
        blacklist=config["arriba"]["blacklist"],
    output:
        fusions="results/arriba/{sample}.fusions.tsv",
        discarded="results/arriva/{sample}.fusions.discarded.tsv",
    threads: 1
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp"
    params:
        blacklist=config["arriba"]["blacklist"],
        extra=config["arriba"].get("extra", ""),
    log:
        "logs/arriba/{sample}.log",
