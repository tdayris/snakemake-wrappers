# Search for fusions with arriba
rule arriba:
    input:
        bam="010.star/{sample}/chimera/{sample}.bam",
        bam_index="010.star/{sample}/chimera/{sample}.bam.bai",
        genome=config["reference"]["genome"],
        genome_index=config["reference"]["genome_index"],
        annotation=config["reference"]["gtf"],
        blacklist=config["arriba"]["blacklist"],
    output:
        fusions=protected("data_output/arriba/{sample}.fusions.tsv"),
        discarded=temp("014.arriba/{sample}.fusions.discarded.tsv"),
    threads: 1
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        blacklist=config["arriba"]["blacklist"],
        extra=config["arriba"].get("extra", ""),
    log:
        "logs/014.arriba/{sample}.log",
    wrapper:
        "bio/arriba"
