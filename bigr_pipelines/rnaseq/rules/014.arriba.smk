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
    params:
        blacklist=config["arriba"]["blacklist"],
        extra="",
    log:
        "logs/arriba/{sample}.log",
