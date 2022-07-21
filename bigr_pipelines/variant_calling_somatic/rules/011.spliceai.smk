rule splice_ai:
    input:
        vcf="bigr/occurence_annotated/{sample}.vcf",
        fasta=config["ref"]["fasta"],
    output:
        vcf=temp("splice_ai/annot/{sample}.vcf.gz"),
    message:
        "Adding Splice Variant annotation to {wildcards.sample}"
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048 * 8,
        time_min=lambda wildcards, attempt: attempt * 60 * 3,
        tmpdir="tmp",
    params:
        annotation=config["params"].get("ncbi_build", "grch38").lower(),
        piped=True,
    log:
        "logs/splice_ai/{sample}.log",
    wrapper:
        "bio/spliceai"
